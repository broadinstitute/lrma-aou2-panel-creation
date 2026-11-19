#!/usr/bin/env bash
# One Google Batch task of the cohort INFO run (direct Batch, no Cromwell). submit.py sets STAGE, CHROM, OUT, RUN, N_BATCHES,
# N_REGIONS; Batch sets BATCH_TASK_INDEX.
#   gp      one task per batch: AF and INFO of every record from the delivered popped BCF's GP (cohort_info.py gp)
#   pool    one task: pool the 200 batches, compare with the Sept re-pop (243ccdd) where it exists, write the cohort INFO
#           onto batch-001 (the paste base)
#   paste   one task per region: slice the base and batches 002-200 by region, paste-vcfs them
#   concat  one task: concatenate the regions and check the merged BCF
# Inputs are only read. Every output goes under $OUT/$CHROM with --no-clobber; a task whose .done exists does nothing.
set -euo pipefail
DELIV=gs://longreadsphase2imputation/deliverables
REPOP=gs://longreadsphase2imputation/deliverables_repop243
O=$OUT/$CHROM
W=/tmp/w; mkdir -p $W; cd $W
gcloud storage cp -q $RUN/cohort_info.py $RUN/paste-vcfs .
chmod +x paste-vcfs

batch() { printf 'batch-%03d' "$1"; }
popped() { echo "$DELIV/$1/aou2_50k_init.$1.$CHROM.prod-05fea60.glimpse2.popped.bcf"; }
put() { gcloud storage cp -q --no-clobber "$@"; }
done_or_exit() { if gcloud storage ls "$1" >/dev/null 2>&1; then echo "already done: $1"; exit 0; fi; }
# region k (1-based) of n equal pieces of the chromosome, from the header of $1
region() {
    local len; len=$(bcftools view -h "$1" | grep "^##contig=<ID=$CHROM," | grep -o 'length=[0-9]*' | cut -d= -f2)
    local step=$(( len / $3 + 1 ))
    echo "$CHROM:$(( ($2 - 1) * step + 1 ))-$(( $2 * step < len ? $2 * step : len ))"
}
INFO_TAGS=ID,RAF,AF,INFO,N_PATHS,N_PATHS_TOTAL,DS_SUM,VAR_SUM,N_HAPS
SITES_INFO_Q='%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/ID\t%INFO/RAF\t%INFO/AF\t%INFO/INFO\t%INFO/N_PATHS\t%INFO/N_PATHS_TOTAL\t%INFO/DS_SUM\t%INFO/VAR_SUM\t%INFO/N_HAPS\n'
GT_Q='%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%GT:%DS:%GP]\n'

stage_gp() {
    local b; b=$(batch $((BATCH_TASK_INDEX + 1)))
    done_or_exit $O/gp/$b.done
    local d; d=$(popped $b)
    gcloud storage cp -q $d $d.csi .
    local f; f=$(basename $d)
    local n; n=$(bcftools query -l $f | awk 'END{print NR}')
    [ "$n" = 250 ] || { echo "ERROR: $b has $n samples" >&2; exit 1; }
    gp_part() {
        set -euo pipefail
        bcftools query -r "$2" --regions-overlap 0 -f '%INFO/INFO[\t%GP]\n' "$3" | python3 cohort_info.py gp > part.$1
    }
    export -f gp_part
    for k in $(seq 1 32); do printf '%02d\t%s\t%s\n' $k "$(region $f $k 32)" $f; done > parts.tsv
    xargs -P "$(nproc)" -L 1 bash -c 'gp_part "$@"' _ < parts.tsv
    cut -f1 parts.tsv | sed 's/^/part./' | xargs cat > info.tsv
    local r; r=$(bcftools index -n $f)
    [ "$(awk 'END{print NR}' info.tsv)" = "$r" ] || { echo "ERROR: $b info rows differ from $r records" >&2; exit 1; }
    gzip -1 info.tsv
    local m; m=$(bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' $f | md5sum | cut -d' ' -f1)
    put info.tsv.gz $O/gp/$b.info.tsv.gz
    printf '%s\t%s\t%s\t%s\n' $b $n $r $m > done
    put done $O/gp/$b.done
}

stage_pool() {
    done_or_exit $O/pool.done
    export GCS_OAUTH_TOKEN; GCS_OAUTH_TOKEN=$(gcloud auth print-access-token)
    local i b
    for i in $(seq 1 $N_BATCHES); do echo $O/gp/$(batch $i).done; done > dones.txt
    gcloud storage cat $(cat dones.txt) > batches.tsv
    [ "$(awk 'END{print NR}' batches.tsv)" = $N_BATCHES ] || { echo "ERROR: not every batch has its gp output" >&2; exit 1; }
    for i in $(seq 1 $N_BATCHES); do echo $O/gp/$(batch $i).info.tsv.gz; done > infos.txt
    mkdir gp; gcloud storage cp -q $(cat infos.txt) gp/

    local r1; r1=$(popped batch-001 | sed "s#^$DELIV#$REPOP#")
    gcloud storage cp -q $r1 $r1.csi .
    bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/ID\t%INFO/RAF\t%INFO/N_PATHS\t%INFO/N_PATHS_TOTAL\n' $(basename $r1) > sites.tsv
    local m; m=$(cut -f1-5 sites.tsv | md5sum | cut -d' ' -f1)
    [ "$(cut -f4 batches.tsv | sort -u)" = "$m" ] || { echo "ERROR: batch site lists differ from each other or from the re-pop" >&2; exit 1; }
    gzip -1 sites.tsv
    python3 cohort_info.py pool --n-samples "$(cut -f2 batches.tsv | paste -sd,)" --sites sites.tsv.gz --report report.txt \
        --out cohort_info.tsv.gz $(cut -f1 batches.tsv | sed 's#.*#gp/&.info.tsv.gz#')

    # the same pool over the batches the Sept re-pop finished, from the GP and from the re-pop's AF/INFO
    : > compare_repop243.txt
    local rb=() rp
    for rp in $(gcloud storage ls "$REPOP/batch-*/*.$CHROM.prod-05fea60.glimpse2.popped.bcf"); do
        b=$(basename "$(dirname $rp)"); rb+=($b)
        bcftools query -f '%INFO/AF\t%INFO/INFO\t.\n' $rp | gzip -1 > repop.$b.info.tsv.gz
    done
    if [ ${#rb[@]} -gt 0 ]; then
        local ns; ns=$(printf '250\n%.0s' "${rb[@]}" | paste -sd,)
        python3 cohort_info.py pool --n-samples $ns --sites sites.tsv.gz --report /dev/null --out sub.gp.tsv.gz \
            $(printf 'gp/%s.info.tsv.gz ' "${rb[@]}")
        python3 cohort_info.py pool --n-samples $ns --sites sites.tsv.gz --report /dev/null --out sub.repop.tsv.gz \
            $(printf 'repop.%s.info.tsv.gz ' "${rb[@]}")
        { echo "batches ${rb[*]}"; echo "a = re-pop 243ccdd pooled, b = GP pooled"
          python3 cohort_info.py compare sub.repop.tsv.gz sub.gp.tsv.gz; } > compare_repop243.txt
    fi

    local d; d=$(popped batch-001)
    gcloud storage cp -q $d $d.csi .
    local base=aou2_50k.cohort_info.$CHROM.base.bcf
    bcftools view $(basename $d) | python3 cohort_info.py apply --table cohort_info.tsv.gz | bcftools view -Ob -o $base --write-index
    [ "$(bcftools index -n $base)" = "$(bcftools index -n $(basename $d))" ] || { echo "ERROR: base record count" >&2; exit 1; }
    [ "$(bcftools query -f "$GT_Q" $base | md5sum)" = "$(bcftools query -f "$GT_Q" $(basename $d) | md5sum)" ] \
        || { echo "ERROR: base sites or genotypes changed" >&2; exit 1; }
    bcftools query -f "$SITES_INFO_Q" $base | md5sum | cut -d' ' -f1 > sites_info.md5
    bcftools index -n $base > n_records
    cat report.txt compare_repop243.txt
    put $base $base.csi cohort_info.tsv.gz report.txt compare_repop243.txt sites_info.md5 n_records batches.tsv $O/pool/
    echo ok > done; put done $O/pool.done
}

stage_paste() {
    local k=$((BATCH_TASK_INDEX + 1)) i
    local tag; tag=$(printf 'region-%03d' $k)
    done_or_exit $O/paste/$tag.done
    gcloud storage ls $O/pool.done >/dev/null
    export GCS_OAUTH_TOKEN; GCS_OAUTH_TOKEN=$(gcloud auth print-access-token)
    local base=$O/pool/aou2_50k.cohort_info.$CHROM.base.bcf
    local reg; reg=$(region "$base##idx##$base.csi" $k $N_REGIONS)
    { echo "001 $base"; for i in $(seq 2 $N_BATCHES); do printf '%03d %s\n' $i "$(popped $(batch $i))"; done; } > inputs.txt
    slice() {
        set -euo pipefail
        local t
        for t in 1 2 3; do
            timeout 30m bcftools view --regions-overlap 0 -r "$REG" "$2##idx##$2.csi" -Ob -o s/$1.bcf --write-index && return 0
            echo "WARNING: retry $t $2" >&2; sleep 10
        done
        exit 1
    }
    export -f slice; export REG=$reg
    mkdir s
    xargs -P 16 -L 1 bash -c 'slice "$@"' _ < inputs.txt
    local want; want=$(bcftools index -n s/001.bcf)
    for f in s/*.bcf; do [ "$(bcftools index -n $f)" = "$want" ] || { echo "ERROR: $f records differ from $want" >&2; exit 1; }; done
    ./paste-vcfs -r $reg --threads "$(nproc)" --info $INFO_TAGS --format GT,DS,GP -o $tag.bcf $(cut -d' ' -f1 inputs.txt | sed 's#.*#s/&.bcf#')
    bcftools index $tag.bcf
    [ "$(bcftools index -n $tag.bcf)" = "$want" ] || { echo "ERROR: pasted records differ from $want" >&2; exit 1; }
    put $tag.bcf $tag.bcf.csi $O/paste/
    printf '%s\t%s\t%s\n' $tag $reg $want > done; put done $O/paste/$tag.done
}

stage_concat() {
    local out=aou2_50k.cohort_info.$CHROM.bcf k
    done_or_exit $O/$out.done
    for k in $(seq 1 $N_REGIONS); do printf '%s/paste/region-%03d.done\n' $O $k; done > dones.txt
    gcloud storage cat $(cat dones.txt) > regions.tsv
    [ "$(awk 'END{print NR}' regions.tsv)" = $N_REGIONS ] || { echo "ERROR: not every region is pasted" >&2; exit 1; }
    gcloud storage cp -q $O/pool/n_records $O/pool/sites_info.md5 $O/pool/batches.tsv .
    mkdir r; gcloud storage cp -q $(cut -f1 regions.tsv | sed "s#.*#$O/paste/&.bcf#") r/
    cut -f1 regions.tsv | sed 's#.*#r/&.bcf#' > list.txt
    bcftools concat --naive -f list.txt -Ob -o $out
    bcftools index --threads "$(nproc)" $out
    rm -r r

    local want_s; want_s=$(cut -f2 batches.tsv | awk '{s += $1} END{print s}')
    local n_s; n_s=$(bcftools query -l $out | awk 'END{print NR}')
    local n_u; n_u=$(bcftools query -l $out | sort -u | awk 'END{print NR}')
    [ "$n_s" = "$want_s" ] && [ "$n_u" = "$n_s" ] || { echo "ERROR: $n_s samples ($n_u unique), expected $want_s" >&2; exit 1; }
    local n_r; n_r=$(bcftools index -n $out)
    [ "$n_r" = "$(cat n_records)" ] || { echo "ERROR: $n_r records, expected $(cat n_records)" >&2; exit 1; }
    [ "$(bcftools query -f "$SITES_INFO_Q" $out | md5sum | cut -d' ' -f1)" = "$(cat sites_info.md5)" ] \
        || { echo "ERROR: merged sites or INFO differ from the base" >&2; exit 1; }
    # genotypes of a few batches over one region, merged vs delivered
    export GCS_OAUTH_TOKEN; GCS_OAUTH_TOKEN=$(gcloud auth print-access-token)
    local reg b d spot=""
    reg=$(region $out 7 $N_REGIONS)
    for b in batch-002 batch-100 batch-200; do
        d=$(popped $b)
        bcftools query -l "$d##idx##$d.csi" > s.$b
        [ "$(bcftools view -S s.$b -r $reg --regions-overlap 0 $out | bcftools query -f "$GT_Q" | md5sum)" \
            = "$(bcftools view -r $reg --regions-overlap 0 "$d##idx##$d.csi" | bcftools query -f "$GT_Q" | md5sum)" ] \
            || { echo "ERROR: $b genotypes in $reg differ from the delivered file" >&2; exit 1; }
        spot="$spot $b"
    done
    printf '%s samples=%s records=%s sites_and_info_equal_to_base=yes genotypes_equal_to_delivered=%s:%s\n' \
        $CHROM "$n_s" "$n_r" "${spot# }" $reg | tee merged_report.txt
    put $out $out.csi merged_report.txt $O/
    echo ok > done; put done $O/$out.done
}

# one task, after every chromosome: copy the merged BCFs, their reports and the README to $DEST (CHROM is unused)
stage_deliver() {
    local c n=0
    gcloud storage ls "$DEST/**" >/dev/null 2>&1 && { echo "ERROR: $DEST is not empty" >&2; exit 1; }
    for c in $(seq 1 22); do gcloud storage ls $OUT/chr$c/aou2_50k.cohort_info.chr$c.bcf.done >/dev/null; done
    gcloud storage cp -q $RUN/README.md README.md
    put README.md $DEST/README.md
    for c in $(seq 1 22); do
        put $OUT/chr$c/aou2_50k.cohort_info.chr$c.bcf $OUT/chr$c/aou2_50k.cohort_info.chr$c.bcf.csi $DEST/
        put $OUT/chr$c/merged_report.txt $OUT/chr$c/pool/report.txt $OUT/chr$c/pool/compare_repop243.txt $DEST/reports/chr$c/
    done
    # same objects, same checksums
    for c in $(seq 1 22); do
        [ "$(gcloud storage ls -L $OUT/chr$c/aou2_50k.cohort_info.chr$c.bcf | grep -i 'crc32c')" = \
          "$(gcloud storage ls -L $DEST/aou2_50k.cohort_info.chr$c.bcf | grep -i 'crc32c')" ] \
            || { echo "ERROR: chr$c copy checksum differs" >&2; exit 1; }
        n=$((n + 1))
    done
    echo "delivered $n chromosomes to $DEST" | tee deliver_report.txt
    put deliver_report.txt $DEST/reports/
}

stage_$STAGE
