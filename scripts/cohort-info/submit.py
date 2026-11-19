#!/usr/bin/env python3
"""Run one chromosome of the cohort INFO fix on Google Batch (no Cromwell): gp -> pool -> paste -> concat, one job per
stage, each submitted when the one before succeeded. Tasks skip work whose .done marker exists, so a rerun resumes.

python3 submit.py CHROM [STAGE ...]

The run files (task.sh, cohort_info.py, paste-vcfs) go to OUT/_run/<md5 of the three>/, so a changed script never
overwrites the one an earlier job ran.
"""
import hashlib
import json
import os
import subprocess
import sys
import time

PROJECT = 'wb-caring-carrot-494'
OUT = 'gs://longreadsphase2imputation/cohort_info_v1'
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
FILES = {'task.sh': os.path.join(HERE, 'task.sh'),
         'cohort_info.py': os.path.join(REPO, 'resources/cohort-info/cohort_info.py'),
         'paste-vcfs': os.path.expanduser('~/aou-info-fix/bin/paste-vcfs'),   # CI-built (utils 243ccdd), SHA256SUMS
         'README.md': os.path.expanduser('~/aou-info-fix/deliverable/README.md')}
DEST = 'gs://long-reads-phase-2-imputation-vwb/cohort_info_v1'
IMAGE = ('us.gcr.io/broad-dsde-methods/slee/lrma-aou2-panel-creation-rust'
         '@sha256:0f25c4091c49d8eb0c3d8bcdb45e7680a093e0a691e74ccec2d3743758a7d22c')
SA = 'pet-27774816681405f7e370a@wb-caring-carrot-494.iam.gserviceaccount.com'
N_BATCHES = 200
# GRCh38. Paste regions of about 3.2 Mb (chr22: 16 regions, about 20 min each).
LENGTHS = {'chr1': 248956422, 'chr2': 242193529, 'chr3': 198295559, 'chr4': 190214555, 'chr5': 181538259,
           'chr6': 170805979, 'chr7': 159345973, 'chr8': 145138636, 'chr9': 138394717, 'chr10': 133797422,
           'chr11': 135086622, 'chr12': 133275309, 'chr13': 114364328, 'chr14': 107043718, 'chr15': 101991189,
           'chr16': 90338345, 'chr17': 83257441, 'chr18': 80373285, 'chr19': 58617616, 'chr20': 64444167,
           'chr21': 46709983, 'chr22': 50818468}
REGION_BP = 3200000
# stage: (tasks, machine, boot disk GB, max run seconds); None = set per chromosome in stages()
STAGES = {'gp': (N_BATCHES, 'n2d-standard-8', 60, 4 * 3600),
          'pool': (1, 'n2d-standard-8', 200, 12 * 3600),
          'paste': (None, 'n2d-standard-8', 200, 12 * 3600),
          'concat': (1, 'n2d-standard-8', None, 12 * 3600),
          'deliver': (1, 'n2d-standard-8', 60, 12 * 3600)}
CHROM_STAGES = ['gp', 'pool', 'paste', 'concat']     # the default; deliver runs once, after all 22: submit.py chr22 deliver


def stages(chrom):
    """Per chromosome: paste region count, and the concat disk from batch-001's size (chr22: 200 popped BCFs 83 GB,
    merged 52 GB; the concat holds the regions and the merged file)."""
    n_regions = -(-LENGTHS[chrom] // REGION_BP)
    f = (f'gs://longreadsphase2imputation/deliverables/batch-001/aou2_50k_init.batch-001.{chrom}'
         '.prod-05fea60.glimpse2.popped.bcf')
    merged_gb = int(sh('gcloud', 'storage', 'du', '-s', f).split()[0]) * N_BATCHES * 0.63 / 1e9
    st = dict(STAGES)
    st['paste'] = (n_regions,) + STAGES['paste'][1:]
    st['concat'] = STAGES['concat'][:2] + (max(600, int(2.4 * merged_gb) + 100),) + STAGES['concat'][3:]
    return n_regions, st


def sh(*cmd):
    return subprocess.run(cmd, capture_output=True, text=True, check=True).stdout


def upload():
    h = hashlib.md5()
    for p in FILES.values():
        h.update(open(p, 'rb').read())
    run = f'{OUT}/_run/{h.hexdigest()[:10]}'
    for name, p in FILES.items():
        sh('gcloud', 'storage', 'cp', '-q', '--no-clobber', p, f'{run}/{name}')
    return run


def job(stage, chrom, run):
    n_regions, st = stages(chrom)
    tasks, machine, disk, max_s = st[stage]
    env = {'STAGE': stage, 'CHROM': chrom, 'OUT': OUT, 'RUN': run, 'DEST': DEST,
           'N_BATCHES': str(N_BATCHES), 'N_REGIONS': str(n_regions)}
    boot = f'mkdir -p /tmp/w && gcloud storage cp -q {run}/task.sh /tmp/w/task.sh && bash /tmp/w/task.sh'
    return {
        'labels': {'stage': f'cohort-info-{stage}', 'chrom': chrom},
        'taskGroups': [{
            'taskCount': str(tasks),
            'parallelism': str(tasks),
            'taskSpec': {
                'computeResource': {'cpuMilli': '8000', 'memoryMib': '28000'},
                'maxRunDuration': f'{max_s}s',
                'maxRetryCount': 3,
                # retry preemptions and VM failures only; a script error fails the task
                'lifecyclePolicies': [{'action': 'FAIL_TASK', 'actionCondition': {'exitCodes': [1, 2, 126, 127, 141]}}],
                'environment': {'variables': env},
                'runnables': [{'container': {'imageUri': IMAGE, 'entrypoint': '/bin/bash', 'commands': ['-c', boot]}}],
            }}],
        'allocationPolicy': {
            'location': {'allowedLocations': ['regions/us-central1']},
            'serviceAccount': {'email': SA},
            'instances': [{'policy': {'machineType': machine, 'provisioningModel': 'SPOT',
                                      'bootDisk': {'type': 'pd-balanced', 'sizeGb': str(disk)}}}],
            'network': {'networkInterfaces': [{
                'network': f'projects/{PROJECT}/global/networks/network',
                'subnetwork': f'projects/{PROJECT}/regions/us-central1/subnetworks/subnetwork',
                'noExternalIpAddress': True}]}},
        'logsPolicy': {'destination': 'CLOUD_LOGGING'},
    }


def submit_and_wait(stage, chrom, run):
    name = f'cohort-info-{chrom}-{stage}-{time.strftime("%m%d%H%M")}'
    path = f'/tmp/{name}.json'
    json.dump(job(stage, chrom, run), open(path, 'w'), indent=1)
    sh('gcloud', 'batch', 'jobs', 'submit', name, '--location=us-central1', f'--project={PROJECT}', f'--config={path}')
    print(f'{time.strftime("%H:%M:%S")} submitted {name}', flush=True)
    while True:
        time.sleep(60)
        try:
            d = json.loads(sh('gcloud', 'batch', 'jobs', 'describe', name, '--location=us-central1',
                              f'--project={PROJECT}', '--format=json'))
        except subprocess.CalledProcessError as e:      # transient API errors; the job keeps running
            print(f'{time.strftime("%H:%M:%S")} describe failed, retrying: {e.stderr.strip()[-200:]}', flush=True)
            continue
        st = d.get('status', {})
        counts = st.get('taskGroups', {}).get('group0', {}).get('counts', {})
        print(f'{time.strftime("%H:%M:%S")} {name} {st.get("state")} {counts}', flush=True)
        if st.get('state') == 'SUCCEEDED':
            return
        if st.get('state') == 'FAILED':
            sys.exit(f'FAILED {name}: see Cloud Logging for job {d.get("uid")}')


def main():
    chrom, stages = sys.argv[1], sys.argv[2:] or CHROM_STAGES
    run = upload()
    print(f'run files {run}', flush=True)
    for s in stages:
        submit_and_wait(s, chrom, run)
    print(f'done {chrom}: {OUT}/{chrom}/aou2_50k.cohort_info.{chrom}.bcf', flush=True)


if __name__ == '__main__':
    main()
