#!/usr/bin/env bash
# check-container-helpers.sh: every helper binary (bcftools, bgzip, tabix,
# samtools) that a script or Nextflow module calls must exist in the image the
# call runs in.
#
# The container smoke test only runs each tool's --version, so it cannot see
# a module that calls bcftools inside the delly image, or a script that runs
# bgzip inside staphb/bcftools, which ships only bcftools.
#
# How the table is built (no hand-kept list):
#   scripts/*.sh   every `docker run ... IMAGE cmd` command, and every call of
#                  the wrapper from scripts/lib/common.sh,
#                  `run_in [--net] [--root] [--rw DIR] ... IMAGE cmd`; the
#                  in-image commands are `cmd`, or, for `bash -c '...'`, the
#                  first word of every command in the quoted body. A pipe
#                  after the closing quote runs on the host and is not counted.
#                  A tree whose scripts yield no helper call at all fails: the
#                  parser no longer sees how the scripts start containers.
#   modules/local/*/main.nf
#                  every process: its `container` (or its withName selector
#                  in conf/containers.config) and the first word of every
#                  command in its script: block.
# Only calls to the helper binaries above are kept. Images are named by their
# versions.env variable when one matches.
#
# Three rules:
#   1. (static) bgzip and tabix are never called in BCFTOOLS_IMAGE or
#      SAMTOOLS_IMAGE: use `bcftools view -Oz -o` and `bcftools index -t`.
#   2. (docker) each image is pulled and `command -v` checks each binary.
#   3. (--orad) orad, Illumina's ORA decompressor, runs on the host, not in an
#      image. The pinned orad release is downloaded, checked against its
#      sha256 and run with --help; every option a script passes to "$ORAD"
#      must be in that list. Step 01 once passed --output-directory, which
#      orad does not have.
#
# Scope: only these four helpers are checked. Probing every in-image command
# (python3, awk, Rscript...) would mean pulling every image the pipeline
# uses, tens of GB, on each pull request.
#
# Usage:
#   scripts/ci/check-container-helpers.sh            both rules (needs docker)
#   scripts/ci/check-container-helpers.sh --static   rule 1 only
#   scripts/ci/check-container-helpers.sh --list     print the derived table
#   scripts/ci/check-container-helpers.sh --orad     rule 3 only (needs curl;
#                                                    --help-file FILE reads
#                                                    the option list from FILE)
#   scripts/ci/check-container-helpers.sh --root DIR check another tree
#   scripts/ci/check-container-helpers.sh --self-test  prove rules 1 and 3
#                                                      can fail, and that a
#                                                      tree whose calls are
#                                                      not recognised fails too
set -euo pipefail

ROOT=$(cd "$(dirname "$0")/../.." && pwd)

if [ "${1:-}" = "--self-test" ]; then
  tmp=$(mktemp -d)
  trap 'rm -rf "$tmp"' EXIT
  for t in bad good; do
    mkdir -p "${tmp}/${t}/scripts" "${tmp}/${t}/modules/local/x"
    printf '%s\n' 'BCFTOOLS_IMAGE="example/bcftools:1.0"' 'SAMTOOLS_IMAGE="example/samtools:1.0"' \
      'TOOL_IMAGE="example/tool:1.0"' > "${tmp}/${t}/versions.env"
  done
  # Planted calls, each behind a different prefix the parser must look past.
  # shellcheck disable=SC2016  # the dollar signs are the test input
  printf '%s\n' \
    '#!/usr/bin/env bash' \
    'if ! docker run --rm -v "$G:/genome" "${BCFTOOLS_IMAGE}" bash -c "bcftools view x | bgzip -c > y.gz && tabix y.gz"; then' \
    '  exit 1' \
    'fi' \
    'OUT=$(time docker run --rm "${BCFTOOLS_IMAGE}" tabix z.gz)' \
    'VAR=x docker run --rm "${SAMTOOLS_IMAGE}" bgzip w' \
    'docker run --rm "${TOOL_IMAGE}" bgzip allowed.vcf' \
    'docker run --rm "${BCFTOOLS_IMAGE}" bash -euo pipefail -c "bgzip combined.vcf"' \
    'IMG="${BCFTOOLS_IMAGE}"' \
    'docker run --rm "$IMG" tabix indirect.vcf.gz' \
    'ODD=$(pick_image)' \
    'docker run --rm "$ODD" samtools view x.bam' \
    'run_in --net --root --rw "$G/reference" --cpus 2 "${SAMTOOLS_IMAGE}" tabix wrapped.vcf.gz' \
    'COUNT=$(run_in "${BCFTOOLS_IMAGE}" bash -c "bcftools view x | bgzip -c | wc -l")' \
    > "${tmp}/bad/scripts/a.sh"
  # A withName selector overrides the process container directive.
  mkdir -p "${tmp}/bad/conf" "${tmp}/bad/modules/local/y"
  printf '%s\n' "process { withName: 'Y' { container = 'example/bcftools:1.0' } }" \
    > "${tmp}/bad/conf/containers.config"
  printf '%s\n' 'process Y {' "    container 'example/tool:1.0'" '    script:' '    """' \
    '    tabix y.vcf.gz' '    """' '}' > "${tmp}/bad/modules/local/y/main.nf"
  printf '%s\n' \
    'process X {' \
    "    container 'example/bcftools:1.0'" \
    '    script:' \
    '    """' \
    '    bcftools view in.vcf | bgzip -c > out.vcf.gz' \
    '    """' \
    '    stub:' \
    '    """' \
    '    tabix stub_only.vcf.gz' \
    '    """' \
    '}' > "${tmp}/bad/modules/local/x/main.nf"
  # shellcheck disable=SC2016
  printf '%s\n' \
    '#!/usr/bin/env bash' \
    'if ! docker run --rm "${BCFTOOLS_IMAGE}" bash -c "bcftools view x -Oz -o y.gz && bcftools index -t y.gz"; then' \
    '  exit 1' \
    'fi' \
    'docker run --rm "${TOOL_IMAGE}" bgzip allowed.vcf' \
    'docker run --rm "${BCFTOOLS_IMAGE}" bash -euo pipefail -c "bcftools index -t ok.vcf.gz"' \
    'run_in --rw "$G/reference" --cpus 2 "${TOOL_IMAGE}" tabix wrapped.vcf.gz' \
    'run_in "${BCFTOOLS_IMAGE}" bcftools index -t wrapped.vcf.gz' \
    > "${tmp}/good/scripts/a.sh"
  # A tree whose scripts start containers in a way the parser does not know:
  # the module alone must not make the check pass.
  mkdir -p "${tmp}/blind/scripts" "${tmp}/blind/modules/local/x"
  cp "${tmp}/good/versions.env" "${tmp}/blind/versions.env"
  # shellcheck disable=SC2016
  printf '%s\n' '#!/usr/bin/env bash' 'start_container "${BCFTOOLS_IMAGE}" bcftools view x' \
    > "${tmp}/blind/scripts/a.sh"
  sed 's/| bgzip -c > out.vcf.gz/-Oz -o out.vcf.gz/' "${tmp}/bad/modules/local/x/main.nf" \
    > "${tmp}/good/modules/local/x/main.nf"

  fail=0
  rc=0; out=$("$0" --static --root "${tmp}/bad" 2>&1) || rc=$?
  [ "$rc" -eq 1 ] || { echo "self-test: planted tree exited ${rc}, expected 1"; fail=1; }
  # shellcheck disable=SC2016  # regexes, not expansions
  for want in 'bgzip +in BCFTOOLS_IMAGE +scripts/a\.sh:2$' 'tabix +in BCFTOOLS_IMAGE +scripts/a\.sh:2$' \
              'tabix +in BCFTOOLS_IMAGE +scripts/a\.sh:5$' 'bgzip +in SAMTOOLS_IMAGE +scripts/a\.sh:6$' \
              'bgzip +in BCFTOOLS_IMAGE +modules/local/x/main\.nf:5 \(X\)$' \
              'bgzip +in BCFTOOLS_IMAGE +scripts/a\.sh:8$' 'tabix +in BCFTOOLS_IMAGE +scripts/a\.sh:10$' \
              'scripts/a\.sh:12 \(\$ODD\)' 'tabix +in BCFTOOLS_IMAGE +modules/local/y/main\.nf:5 \(Y\)$' \
              'tabix +in SAMTOOLS_IMAGE +scripts/a\.sh:13$' 'bgzip +in BCFTOOLS_IMAGE +scripts/a\.sh:14$'; do
    grep -qE -- "$want" <<<"$out" || { echo "self-test: planted call not reported: ${want}"; fail=1; }
  done
  for nope in 'TOOL_IMAGE' 'stub_only' ':9 \(X\)'; do
    if grep -qE -- "$nope" <<<"$out"; then echo "self-test: reported although allowed: ${nope}"; fail=1; fi
  done
  [ "$fail" -eq 0 ] || printf '%s\n' "$out"
  rc=0; out=$("$0" --static --root "${tmp}/good" 2>&1) || rc=$?
  if [ "$rc" -ne 0 ]; then
    echo "self-test: clean tree exited ${rc}, expected 0"
    printf '%s\n' "$out"
    fail=1
  fi
  sed 's/| bgzip -c > out.vcf.gz/-Oz -o out.vcf.gz/' "${tmp}/bad/modules/local/x/main.nf" \
    > "${tmp}/blind/modules/local/x/main.nf"
  rc=0; out=$("$0" --static --root "${tmp}/blind" 2>&1) || rc=$?
  if [ "$rc" -ne 2 ] || ! grep -q 'no helper call in scripts/' <<<"$out"; then
    echo "self-test: a tree whose script calls are not recognised exited ${rc}, expected 2 with 'no helper call in scripts/'"
    printf '%s\n' "$out"
    fail=1
  fi
  # Rule 3 against a planted orad option list: an option orad does not list
  # fails, the options it lists pass, and a tree with no orad call or an
  # empty help text fails too.
  printf '%s\n' '-P  path, --path  path   : specify where to write the output file' \
    '-h, --help               : print help and exit' \
    '--ora-reference          : set directory path containing the ora reference file' > "${tmp}/orad-help"
  : > "${tmp}/orad-help-empty"
  # shellcheck disable=SC2016,SC1003  # literal $ and line-continuation backslashes are the test input
  printf '%s\n' '#!/usr/bin/env bash' 'if ! command -v "$ORAD" >/dev/null; then exit 1; fi' \
    '"$ORAD" \' '  --ora-reference "$R" \' '  --output-directory "$O" \' '  "$F"' > "${tmp}/bad/scripts/b.sh"
  # shellcheck disable=SC2016
  printf '%s\n' '#!/usr/bin/env bash' 'time "${ORAD}" --ora-reference="$R" -P "$O" --path "$O" "$F"' \
    > "${tmp}/good/scripts/b.sh"
  rc=0; out=$("$0" --orad --help-file "${tmp}/orad-help" --root "${tmp}/bad" 2>&1) || rc=$?
  if [ "$rc" -ne 1 ] || ! grep -q 'scripts/b\.sh:3 passes --output-directory' <<<"$out" \
      || ! grep -q 'OK: --ora-reference' <<<"$out"; then
    echo "self-test: an orad option missing from --help exited ${rc}, expected 1 naming scripts/b.sh:3 --output-directory"
    printf '%s\n' "$out"
    fail=1
  fi
  rc=0; out=$("$0" --orad --help-file "${tmp}/orad-help" --root "${tmp}/good" 2>&1) || rc=$?
  if [ "$rc" -ne 0 ] || [ "$(grep -c '^OK: ' <<<"$out")" -ne 3 ]; then
    echo "self-test: orad options that --help lists exited ${rc}, expected 0 with 3 OK lines"
    printf '%s\n' "$out"
    fail=1
  fi
  rc=0; out=$("$0" --orad --help-file "${tmp}/orad-help-empty" --root "${tmp}/good" 2>&1) || rc=$?
  if [ "$rc" -ne 1 ]; then
    echo "self-test: an empty orad --help exited ${rc}, expected 1"
    printf '%s\n' "$out"
    fail=1
  fi
  rc=0; out=$("$0" --orad --help-file "${tmp}/orad-help" --root "${tmp}/blind" 2>&1) || rc=$?
  # shellcheck disable=SC2016  # a literal $ORAD in the expected message
  if [ "$rc" -ne 2 ] || ! grep -q 'no "\$ORAD" call' <<<"$out"; then
    echo "self-test: a tree with no orad call exited ${rc}, expected 2"
    printf '%s\n' "$out"
    fail=1
  fi
  if [ "$fail" -eq 0 ]; then
    echo "self-test: bgzip/tabix behind if !, time, VAR=, \$(...), bash -euo -c, an image alias, the run_in wrapper and in modules are reported; a clean tree passes; a tree whose script calls are not recognised fails; an orad option missing from --help, an empty --help and a tree with no orad call fail: PASS"
    exit 0
  fi
  echo "self-test: FAIL"
  exit 1
fi

exec python3 - "$ROOT" "$@" <<'PY'
import glob
import os
import re
import subprocess
import sys

ROOT = sys.argv[1]
ARGS = sys.argv[2:]
if "--root" in ARGS:
    ROOT = os.path.abspath(ARGS[ARGS.index("--root") + 1])
HELPERS = ("bcftools", "bgzip", "tabix", "samtools")
NO_BGZIP_TABIX = ("BCFTOOLS_IMAGE", "SAMTOOLS_IMAGE")
KEYWORDS = {"if", "then", "elif", "else", "fi", "do", "done", "while", "until",
            "{", "}", "!", "time", "function", "exec", "[[", "]]"}
STOP = {"for", "select", "case", "in", "esac"}

# ---------------------------------------------------------------- versions
def versions():
    script = ('set -euo pipefail; . "$1"; for v in $(compgen -v); do '
              'case $v in *_IMAGE) printf "%s=%s\\n" "$v" "${!v}";; esac; done')
    out = subprocess.run(["env", "-i", "bash", "--noprofile", "--norc", "-c", script,
                          "_", os.path.join(ROOT, "versions.env")],
                         capture_output=True, text=True, check=True).stdout
    return dict(line.split("=", 1) for line in out.splitlines() if "=" in line)

VERSIONS = versions()
BY_VALUE = {v: k for k, v in VERSIONS.items()}

# ---------------------------------------------------------------- shell lexer
def skip_dq(t, j):
    """t[j] is a double quote; return the index after its closing quote."""
    j += 1
    while j < len(t) and t[j] != '"':
        if t[j] == "\\":
            j += 2
            continue
        if t.startswith("$(", j):
            j = match_paren(t, j + 1) + 1
            continue
        j += 1
    return j + 1

def match_paren(t, i):
    """t[i] is '('; return the index of its matching ')'."""
    depth, j = 0, i
    while j < len(t):
        c = t[j]
        if c == "\\":
            j += 2
            continue
        if c == "'":
            k = t.find("'", j + 1)
            j = len(t) if k < 0 else k + 1
            continue
        if c == '"':
            j = skip_dq(t, j)
            continue
        if c == "(":
            depth += 1
        elif c == ")":
            depth -= 1
            if depth == 0:
                return j
        j += 1
    return len(t) - 1

def match_brace(t, i):
    """t[i] is '{'; return the index of its matching '}'."""
    depth, j = 0, i
    while j < len(t):
        if t[j] == "{":
            depth += 1
        elif t[j] == "}":
            depth -= 1
            if depth == 0:
                return j
        j += 1
    return len(t) - 1

def lex(text, base=0):
    """Split shell text into simple commands: (words, offset of first word).
    Command and process substitutions are lexed too and their commands added;
    here-document bodies, comments and redirection targets are dropped."""
    cmds, cur, heredocs = [], [], []
    word, wstart, redirect = None, 0, False
    i, n = 0, len(text)

    def add(s, at):
        nonlocal word, wstart
        if word is None:
            word, wstart = "", at
        word += s

    def end_word():
        nonlocal word, redirect
        if word is not None:
            if redirect:
                redirect = False
            else:
                cur.append((word, wstart))
        word = None

    def end_cmd():
        nonlocal cur
        end_word()
        if cur:
            cmds.append(([w for w, _ in cur], base + cur[0][1]))
        cur = []

    while i < n:
        c = text[i]
        if c == "\\":
            if i + 1 < n and text[i + 1] == "\n":
                i += 2
                continue
            add(text[i + 1:i + 2], i)
            i += 2
        elif c == "'":
            k = text.find("'", i + 1)
            k = n if k < 0 else k
            add(text[i + 1:k], i)
            i = k + 1
        elif c == '"':
            j, buf = i + 1, ""
            while j < n and text[j] != '"':
                if text[j] == "\\" and j + 1 < n and text[j + 1] in '"$\\`\n':
                    if text[j + 1] != "\n":
                        buf += text[j + 1]
                    j += 2
                    continue
                if text.startswith("$(", j) and not text.startswith("$((", j):
                    k = match_paren(text, j + 1)
                    cmds.extend(lex(text[j + 2:k], base + j + 2))
                    buf += "$(...)"
                    j = k + 1
                    continue
                buf += text[j]
                j += 1
            add(buf, i)
            i = j + 1
        elif (c == "$" and text.startswith("$(", i) and not text.startswith("$((", i)) or \
                (c in "<>" and text.startswith("(", i + 1)):
            k = match_paren(text, i + 1)
            cmds.extend(lex(text[i + 2:k], base + i + 2))
            add("$(...)", i)
            i = k + 1
        elif c == "$" and text.startswith("${", i):
            k = match_brace(text, i + 1)
            add(text[i:k + 1], i)
            i = k + 1
        elif c == "#" and word is None:
            k = text.find("\n", i)
            i = n if k < 0 else k
        elif c in " \t":
            end_word()
            i += 1
        elif c == "\n":
            end_cmd()
            i += 1
            for delim, strip in heredocs:
                while i < n:
                    k = text.find("\n", i)
                    k = n if k < 0 else k
                    line = text[i:k]
                    i = k + 1
                    if (line.lstrip("\t") if strip else line).strip() == delim:
                        break
            heredocs = []
        elif c in ";|&()":
            if text.startswith("&>", i):
                end_word()
                redirect = True
                i += 2
            elif text[i:i + 2] in ("||", "&&", ";;", "|&"):
                end_cmd()
                i += 2
            else:
                end_cmd()
                i += 1
        elif c in "<>":
            if text.startswith("<<<", i):
                end_word()
                i += 3
            elif text.startswith("<<", i):
                end_word()
                j, strip = i + 2, False
                if j < n and text[j] == "-":
                    strip, j = True, j + 1
                while j < n and text[j] in " \t":
                    j += 1
                m = re.match(r"""(['"]?)([A-Za-z0-9_]+)\1""", text[j:])
                if m:
                    heredocs.append((m.group(2), strip))
                    j += m.end()
                i = j
            else:
                if word is not None and word.isdigit():
                    word = None
                end_word()
                j = i
                while j < n and text[j] in "<>&|":
                    j += 1
                redirect = True
                i = j
        else:
            add(c, i)
            i += 1
    end_cmd()
    return cmds

ASSIGN = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*(\[[^]]*\])?\+?=")

def heads(cmds):
    """(binary, words, offset) for the command word of every simple command."""
    out = []
    for words, off in cmds:
        k = 0
        while k < len(words):
            w = words[k]
            if w in STOP:
                break
            if w in KEYWORDS or ASSIGN.match(w):
                k += 1
                continue
            out.append((os.path.basename(w), words[k:], off))
            break
    return out

# ---------------------------------------------------------------- docker run
VALUED = set("""-a --attach --add-host -c --cpu-shares --cap-add --cap-drop --cidfile --cpus
 --cpuset-cpus --device --dns -e --env --env-file --entrypoint --expose --gpus --group-add
 -h --hostname --ipc -l --label --label-file --link --log-driver --log-opt -m --memory
 --memory-reservation --memory-swap --mount --name --net --network -p --publish --pid
 --platform --pull --restart --runtime --security-opt --shm-size --stop-signal --tmpfs
 -u --user --ulimit --userns --uts -v --volume --volumes-from -w --workdir""".split())

def docker_run(words):
    """(image word, in-image argv) for `docker run ...` or for the wrapper
    `run_in [--net] [--root] [--rw DIR]... [docker run options] ...`, or None."""
    if len(words) >= 2 and words[0] == "docker" and words[1] == "run":
        k = 2
    elif words and words[0] == "run_in":
        k = 1
        while k < len(words) and words[k] in ("--net", "--root", "--rw"):
            k += 2 if words[k] == "--rw" else 1
    else:
        return None
    entry = []
    while k < len(words) and words[k].startswith("-"):
        w = words[k]
        if w == "--entrypoint" and k + 1 < len(words):
            entry = [words[k + 1]]
        k += 2 if (w in VALUED and "=" not in w) else 1
    if k >= len(words):
        return None
    return words[k], entry + words[k + 1:]

def in_image_heads(argv):
    """Binaries a container command runs: argv[0], or every command of a
    `bash -c` / `sh -c` body. Option groups such as `-euo pipefail` are
    understood: each `o` or `O` takes the next word, and a `c` anywhere in a
    group makes the first word after the options the command body."""
    if not argv:
        return []
    if os.path.basename(argv[0]) not in ("bash", "sh"):
        return [os.path.basename(argv[0])]
    k, want_c, takes = 1, False, 0
    while k < len(argv):
        w = argv[k]
        if takes:
            takes -= 1
        elif w == "--":
            k += 1
            break
        elif w.startswith("--"):
            pass
        elif len(w) > 1 and w[0] in "-+":
            want_c = want_c or "c" in w[1:]
            takes = w[1:].count("o") + w[1:].count("O")
        else:
            break
        k += 1
    if want_c and k < len(argv):
        return [h for h, _, _ in heads(lex(argv[k]))]
    return []

PLAIN = re.compile(r"^[A-Za-z0-9._/:@+-]+$")

def resolve(image_word, path):
    """(variable name, image value) for an image word from a script. A
    script-local assignment counts only when its value is a plain literal or
    one versions.env variable; anything else is unresolved (value None)."""
    m = re.search(r"\$\{?([A-Za-z_][A-Za-z0-9_]*)", image_word)
    if m:
        name = m.group(1)
        if name in VERSIONS:
            return name, VERSIONS[name]
        with open(path) as fh:
            for line in fh:
                a = re.match(r"\s*(?:export\s+|local\s+|readonly\s+)?%s=(\S*)" % re.escape(name), line)
                if not a:
                    continue
                v = a.group(1).rstrip(";")
                if len(v) >= 2 and v[0] == v[-1] and v[0] in "\"'":
                    v = v[1:-1]
                ref = re.match(r"^\$\{?([A-Za-z_][A-Za-z0-9_]*)\}?$", v)
                if ref and ref.group(1) in VERSIONS:
                    return ref.group(1), VERSIONS[ref.group(1)]
                if PLAIN.match(v):
                    return name, v
                return name, None
        return name, None
    return BY_VALUE.get(image_word, image_word), image_word

def line_of(text, offset):
    return text.count("\n", 0, offset) + 1

# ---------------------------------------------------------------- orad
# Rule 3. orad is Illumina's own binary, not an image: step 01 calls it on the
# host. The pinned release is fetched (only its first 9 MB: the binary comes
# before the 785 MB reference in the archive), its checksum is checked, and
# the options of every `"$ORAD" ...` call in scripts/*.sh must appear in its
# --help. --help-file FILE reads the option list from FILE instead.
ORAD_VERSION = "2.7.0"
ORAD_URL = ("https://webdata.illumina.com/downloads/software/dragen-decompression/"
            "orad.%s.linux.tar.gz" % ORAD_VERSION)
ORAD_MEMBER = "orad.%s.linux/orad" % ORAD_VERSION
ORAD_SHA256 = "2cc4a6bb9d0721c6555e01bb23f3554753166bfdbdd89ee9f4a8639144a1f066"

def orad_calls():
    """(location, options) for every `"$ORAD" ...` command in scripts/*.sh."""
    found = []
    for path in sorted(glob.glob(os.path.join(ROOT, "scripts", "*.sh"))):
        with open(path) as fh:
            text = fh.read()
        for words, off in lex(text):
            k = 0
            while k < len(words) and (words[k] in KEYWORDS or ASSIGN.match(words[k])):
                k += 1
            if k < len(words) and re.match(r"^\$\{?ORAD\}?$", words[k]):
                opts = [w.split("=", 1)[0] for w in words[k + 1:] if len(w) > 1 and w[0] == "-"]
                found.append(("%s:%d" % (os.path.relpath(path, ROOT), line_of(text, off)), opts))
    return found

def orad_help():
    """orad's --help text, or None with the reason printed."""
    if "--help-file" in ARGS:
        with open(ARGS[ARGS.index("--help-file") + 1]) as fh:
            return fh.read()
    import hashlib, shlex, tempfile
    tmp = tempfile.mkdtemp()
    gnu = "GNU" in subprocess.run(["tar", "--version"], capture_output=True, text=True).stdout
    stop = "--occurrence=1" if gnu else "-q"   # stop reading once the binary is out
    # The exit status is not used: tar stops early, so curl ends on a broken
    # pipe. The checksum below is what proves the right binary arrived.
    # --retry covers a transient network error; curl never retries its (23)
    # write error, so the early stop still ends the fetch.
    print("Fetching %s from %s (tar stops after it, so a curl (23) write error"
          " below is expected)" % (ORAD_MEMBER, ORAD_URL), flush=True)
    subprocess.run("curl -fsSL --retry 3 %s | tar -xzf - -C %s %s %s" % (
        shlex.quote(ORAD_URL), shlex.quote(tmp), stop, shlex.quote(ORAD_MEMBER)), shell=True)
    binary = os.path.join(tmp, ORAD_MEMBER)
    if not os.path.isfile(binary):
        print("FAIL: could not fetch %s from %s" % (ORAD_MEMBER, ORAD_URL))
        return None
    with open(binary, "rb") as fh:
        sha = hashlib.sha256(fh.read()).hexdigest()
    if sha != ORAD_SHA256:
        print("FAIL: orad %s has sha256 %s, pinned %s" % (ORAD_VERSION, sha, ORAD_SHA256))
        return None
    os.chmod(binary, 0o755)
    res = subprocess.run([binary, "--help"], capture_output=True, text=True)
    return res.stdout + res.stderr

if "--orad" in ARGS:
    calls = orad_calls()
    if not calls:
        print("ERROR: found no \"$ORAD\" call in scripts/*.sh; the parser does not see how step 01 runs orad.")
        sys.exit(2)
    text = orad_help()
    if text is None:
        sys.exit(1)
    # An empty or foreign help text must not pass as "every option listed".
    if not re.search(r"^\s*-h, --help\b", text, re.M):
        print("FAIL: orad --help printed no option list:")
        print(text)
        sys.exit(1)
    failed = False
    for loc, opts in calls:
        for o in opts:
            if re.search(r"(?<![\w-])%s(?![\w-])" % re.escape(o), text):
                print("OK: %-20s %s" % (o, loc))
            else:
                failed = True
                print("FAIL: %s passes %s, which orad %s does not list in --help" % (loc, o, ORAD_VERSION))
    if failed:
        print("orad %s --help:" % ORAD_VERSION)
        print(text)
    sys.exit(1 if failed else 0)

# ---------------------------------------------------------------- collect
calls = []        # (image name, image value, binary, location)
unresolved = []   # locations whose image could not be determined

for path in sorted(glob.glob(os.path.join(ROOT, "scripts", "*.sh"))):
    rel = os.path.relpath(path, ROOT)
    with open(path) as fh:
        text = fh.read()
    for words, off in lex(text):
        # `if ! docker run`, `time docker run`, `VAR=x docker run`: skip the
        # keywords and assignments in front, as heads() does.
        k = 0
        while k < len(words) and (words[k] in KEYWORDS or ASSIGN.match(words[k])):
            k += 1
        dr = docker_run(words[k:])
        if not dr:
            continue
        image_word, argv = dr
        name, value = resolve(image_word, path)
        for b in in_image_heads(argv):
            if b in HELPERS:
                hit = re.compile(r"(?<![\w/-])%s(?![\w-])" % re.escape(b)).search(text, off)
                loc = "%s:%d" % (rel, line_of(text, hit.start() if hit else off))
                if value is None:
                    unresolved.append("%s (%s)" % (loc, image_word))
                else:
                    calls.append((name, value, b, loc))

def groovy_to_shell(s, interpolate):
    out, i = [], 0
    while i < len(s):
        c = s[i]
        if c == "\\" and i + 1 < len(s):
            nxt = s[i + 1]
            out.append(nxt if nxt in "$\\\"'" else c + nxt)
            i += 2
            continue
        if c == "$" and interpolate:
            if s.startswith("${", i):
                i = match_brace(s, i + 1) + 1
                out.append("GROOVY")
                continue
            m = re.match(r"\$[A-Za-z_][A-Za-z0-9_.]*", s[i:])
            if m:
                i += m.end()
                out.append("GROOVY")
                continue
        out.append(c)
        i += 1
    return "".join(out)

selectors = {}
cfg = os.path.join(ROOT, "conf", "containers.config")
if os.path.exists(cfg):
    with open(cfg) as fh:
        for m in re.finditer(r"withName:\s*['\"]?(\w+)['\"]?\s*\{[^}]*?container\s*=\s*['\"]([^'\"]+)['\"]",
                             fh.read()):
            selectors[m.group(1)] = m.group(2)

PROC = re.compile(r"^process\s+(\w+)\s*\{", re.M)
for path in sorted(glob.glob(os.path.join(ROOT, "modules", "local", "*", "main.nf"))):
    rel = os.path.relpath(path, ROOT)
    with open(path) as fh:
        text = fh.read()
    procs = list(PROC.finditer(text))
    for idx, m in enumerate(procs):
        pname, start = m.group(1), m.end()
        end = procs[idx + 1].start() if idx + 1 < len(procs) else len(text)
        body = text[start:end]
        cm = re.search(r"^\s*container\s+['\"]([^'\"]+)['\"]", body, re.M)
        # A withName selector in the config overrides the process directive.
        image = selectors.get(pname) or (cm.group(1) if cm else None)
        sm = re.search(r"^\s*script:\s*$", body, re.M)
        if not sm:
            continue
        sec_start = start + sm.end()
        stub = re.search(r"^\s*stub:\s*$", text[sec_start:end], re.M)
        sec_end = sec_start + stub.start() if stub else end
        for q in re.finditer(r'("""|\'\'\')(.*?)\1', text[sec_start:sec_end], re.S):
            shell = groovy_to_shell(q.group(2), q.group(1) == '"""')
            qstart = sec_start + q.start(2)
            for b, _, _ in heads(lex(shell)):
                if b not in HELPERS:
                    continue
                hit = re.search(r"(?<![\w/-])%s(?![\w-])" % re.escape(b), text[qstart:sec_end])
                loc = "%s:%d (%s)" % (rel, line_of(text, qstart + (hit.start() if hit else 0)), pname)
                if image is None:
                    unresolved.append(loc)
                else:
                    calls.append((BY_VALUE.get(image, image), image, b, loc))

# de-duplicate, keep order
seen, table = set(), []
for c in calls:
    if c not in seen:
        seen.add(c)
        table.append(c)

# The scripts call bcftools and samtools in containers all over. If none of
# those calls was derived, the parser no longer recognises how the scripts
# start a container, and a green result would mean nothing.
has_scripts = bool(glob.glob(os.path.join(ROOT, "scripts", "*.sh")))
blind = has_scripts and not any(loc.startswith("scripts/") for _, _, _, loc in table) \
    and not any(loc.startswith("scripts/") for loc in unresolved)

if "--list" in ARGS:
    for name, value, b, loc in table:
        print("%-16s %-9s %s" % (name, b, loc))
    if blind:
        print("ERROR: derived no helper call in scripts/*.sh; the parser does not see how the scripts start containers.")
        sys.exit(2)
    sys.exit(0)

if blind:
    print("ERROR: derived no helper call in scripts/*.sh; the parser does not see how the scripts start containers.")
    sys.exit(2)

failed = False
if not table:
    print("ERROR: derived no helper calls at all; the parser is broken.")
    sys.exit(2)
if unresolved:
    failed = True
    print("FAIL: helper calls whose image could not be determined:")
    for loc in unresolved:
        print("  " + loc)

bad = [(n, b, loc) for n, _, b, loc in table if n in NO_BGZIP_TABIX and b in ("bgzip", "tabix")]
if bad:
    failed = True
    print("FAIL: bgzip/tabix called in an image that is not meant to provide them")
    print("      (use `bcftools view -Oz -o X.vcf.gz` and `bcftools index -t X.vcf.gz`):")
    for n, b, loc in bad:
        print("  %-9s in %-15s %s" % (b, n, loc))
else:
    print("OK: no bgzip/tabix call in %s." % " or ".join(NO_BGZIP_TABIX))

if "--static" in ARGS:
    sys.exit(1 if failed else 0)

images = {}
for name, value, b, loc in table:
    images.setdefault((name, value), {}).setdefault(b, []).append(loc)

missing_total = 0
for (name, value), bins in sorted(images.items()):
    print("::group::%s (%s): %s" % (name, value, " ".join(sorted(bins))), flush=True)
    pull = subprocess.run(["docker", "pull", "-q", value], capture_output=True, text=True)
    if pull.returncode != 0:
        print("::endgroup::")
        print("FAIL: cannot pull %s: %s" % (value, pull.stderr.strip()))
        failed = True
        continue
    probe = ('for b in "$@"; do command -v "$b" >/dev/null 2>&1 && echo "found $b" '
             '|| echo "missing $b"; done')
    res = subprocess.run(["docker", "run", "--rm", "--entrypoint", "sh", value, "-c", probe, "sh"]
                         + sorted(bins), capture_output=True, text=True)
    print(res.stdout + res.stderr, end="")
    print("::endgroup::")
    found = set(re.findall(r"^found (\S+)$", res.stdout, re.M))
    if res.returncode != 0 or not (found or "missing" in res.stdout):
        print("FAIL: could not probe %s (exit %d)" % (value, res.returncode))
        failed = True
        continue
    for b in sorted(bins):
        if b not in found:
            missing_total += 1
            failed = True
            print("FAIL: %s is not in %s (%s), called at:" % (b, name, value))
            for loc in bins[b]:
                print("    " + loc)
    subprocess.run(["docker", "rmi", "-f", value], capture_output=True)

if missing_total == 0 and not failed:
    print("OK: every helper binary is present in the image that calls it (%d images)." % len(images))
sys.exit(1 if failed else 0)
PY
