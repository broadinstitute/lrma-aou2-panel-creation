#!/usr/bin/env python3
"""
Static checks on the shell inside each WDL task command block.

The argument/sizing tests exercise a *mirror* of the command block, so they cannot see the
real thing. This parses the WDL directly and catches what the mirror cannot: a call to a
shell function that is not defined in that block (under set -e that is exit 127 on every
task), a syntax error, or a [RESOURCE] line written to stdout instead of stderr.

    python3 test/wdl/check-command-blocks.py <file.wdl>
"""
import re, subprocess, sys, tempfile, os
wdl = sys.argv[1]
s = open(wdl).read()
fail = 0

# Commands legitimately invoked in these blocks. Anything else that looks like a shell
# helper (snake_case, or leading underscore) and is not defined in the same command block
# is a call to a function that does not exist -- which under set -e is exit 127.
KNOWN = {
    "set","echo","printf","cat","grep","sed","awk","paste","cut","wc","df","kill","trap",
    "exit","eval","cp","mv","rm","mkdir","chmod","touch","read","while","for","if","then",
    "else","elif","fi","do","done","case","esac","local","return","export","bcftools",
    "cargo","sleep","true","false","tail","head","sort","uniq","mktemp","dirname","basename",
    # binaries provided by the task images, invoked bare rather than by path
    "ligate_static",
}
HELPERISH = re.compile(r"^_[A-Za-z0-9_]+$|^[a-z][a-z0-9]*_[a-z0-9_]+$")

for m in re.finditer(r"task (\w+) \{", s):
    name = m.group(1)
    seg = s[m.end():]
    nx = re.search(r"\ntask ", seg)
    seg = seg[:nx.start()] if nx else seg
    cm = re.search(r"command <<<(.*?)\n\s*>>>", seg, re.S)
    if not cm:
        continue
    body = cm.group(1)

    # Strip embedded heredocs (python/awk/etc.) before looking for shell invocations --
    # a Python line like "bin_df = ..." is not a command call.
    shell_only = re.sub(r"<<-?\s*[\"']?(\w+)[\"']?\n.*?\n\s*\1\s*$", "", body,
                        flags=re.S | re.M)

    local = set(re.findall(r"^\s*([A-Za-z_][A-Za-z0-9_]*)\(\)\s*\{", shell_only, re.M))
    invoked = set()
    for line in shell_only.splitlines():
        t = line.strip()
        if not t or t.startswith("#"):
            continue
        # a command invocation: NAME followed by whitespace and something that is not '=',
        # which excludes both "NAME=value" and "NAME = value" (assignment in any language)
        tok = re.match(r"([A-Za-z_][A-Za-z0-9_]*)\s+(?![=|&])", t)
        if tok and not re.match(r"[A-Za-z_][A-Za-z0-9_]*\s*=", t):
            invoked.add(tok.group(1))
    suspects = sorted(x for x in invoked
                      if x not in KNOWN and x not in local and HELPERISH.match(x))
    if suspects:
        print("  FAIL %s: calls undefined shell function(s): %s" % (name, ", ".join(suspects)))
        fail = 1

    shell = re.sub(r"~\{[^}]*\}", "WDLPLACEHOLDER", body)
    with tempfile.NamedTemporaryFile("w", suffix=".sh", delete=False) as fh:
        fh.write(shell); path = fh.name
    r = subprocess.run(["bash", "-n", path], capture_output=True, text=True)
    os.unlink(path)
    if r.returncode != 0:
        print("  FAIL %s: bash -n: %s" % (name, r.stderr.strip().splitlines()[:2]))
        fail = 1

    for line in body.splitlines():
        if "[RESOURCE]" in line and "echo" in line and not line.rstrip().endswith(">&2"):
            print("  FAIL %s: [RESOURCE] echo not redirected to stderr" % name)
            fail = 1
            break

if not fail:
    print("  all command blocks parse, define what they call, and log to stderr")
sys.exit(fail)
