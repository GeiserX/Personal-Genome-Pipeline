#!/usr/bin/env bash
# scripts/ci/check-exec-bits.sh on a planted tree: a helper a module calls by
# name, committed without its executable bit, fails the check; with the bit
# it passes; a helper named only as a path (bin/x.py in a comment) is not a
# call; an executable helper without a #! line fails. Then the real tree.
set -uo pipefail

REPO="$(cd "$(dirname "$0")/.." && pwd)"
CHECK="${REPO}/scripts/ci/check-exec-bits.sh"
T=$(mktemp -d)
trap 'rm -rf "$T"' EXIT
fails=0
ok()  { echo "PASS $*"; }
bad() { echo "FAIL $*"; fails=$((fails + 1)); }

g() { git -C "$T" -c user.name=t -c user.email=t@example.invalid "$@"; }
g init -q
mkdir -p "${T}/bin" "${T}/modules/local/scan"
printf '#!/usr/bin/env python3\nprint(1)\n' > "${T}/bin/scan_it.py"
printf '#!/usr/bin/env python3\nprint(2)\n' > "${T}/bin/only_named.py"
cat > "${T}/modules/local/scan/main.nf" <<'NF'
process SCAN {
    // reads what bin/only_named.py wrote
    script:
    """
    scan_it.py --out x.tsv
    """
}
NF
printf 'workflow { SCAN() }\n' > "${T}/main.nf"
g add -A
g update-index --chmod=-x bin/scan_it.py bin/only_named.py

out=$(bash "$CHECK" "$T" 2>&1); rc=$?
if [ "$rc" -ne 0 ] && grep -q 'FAIL bin/scan_it.py: mode 100644' <<< "$out"; then ok "a called helper without the bit fails"
else bad "a called helper without the bit fails (exit ${rc}): ${out}"; fi
if grep -q 'only_named' <<< "$out"; then bad "a helper named only as a path is not a call: ${out}"
else ok "a helper named only as a path is not a call"; fi

g update-index --chmod=+x bin/scan_it.py
out=$(bash "$CHECK" "$T" 2>&1); rc=$?
if [ "$rc" -eq 0 ]; then ok "with the bit set it passes"; else bad "with the bit set it passes (exit ${rc}): ${out}"; fi

printf 'print(1)\n' > "${T}/bin/scan_it.py"
out=$(bash "$CHECK" "$T" 2>&1); rc=$?
if [ "$rc" -ne 0 ] && grep -q 'not a #! line' <<< "$out"; then ok "an executable helper without #! fails"
else bad "an executable helper without #! fails (exit ${rc}): ${out}"; fi

printf 'workflow { }\n' > "${T}/main.nf"
rm -rf "${T}/modules"
out=$(bash "$CHECK" "$T" 2>&1); rc=$?
if [ "$rc" -ne 0 ]; then ok "a tree where no helper is called fails (the parser saw nothing)"
else bad "a tree where no helper is called fails: ${out}"; fi

out=$(bash "$CHECK" "$REPO" 2>&1); rc=$?
echo "$out" | tail -n 3
if [ "$rc" -eq 0 ]; then ok "this checkout"; else bad "this checkout (exit ${rc})"; fi

[ "$fails" -eq 0 ] || { echo "${fails} check(s) failed"; exit 1; }
echo "all checks passed"
