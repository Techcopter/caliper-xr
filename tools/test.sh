#!/usr/bin/env bash
# Run the headless test suite against this checkout.
#   tools/test.sh                  all suites
#   tools/test.sh tests/test_shell.py
# BLENDER overrides the Blender binary (default: the macOS app bundle).
set -u
BLENDER="${BLENDER:-/Applications/Blender.app/Contents/MacOS/Blender}"
cd "$(dirname "$0")/.."
suites=("$@")
[ ${#suites[@]} -eq 0 ] && suites=(tests/test_*.py)
failed=()
for t in "${suites[@]}"; do
  printf '\n### %s\n' "$t"
  "$BLENDER" -b --factory-startup --python-exit-code 1 --python "$t" 2>&1 \
    | grep -E '^(  PASS|  FAIL|  info|== )|Error|Traceback|  File "'
  [ "${PIPESTATUS[0]}" -eq 0 ] || failed+=("$t")
done
echo
if [ ${#failed[@]} -eq 0 ]; then
  echo "ALL SUITES PASSED (${#suites[@]})"
else
  echo "FAILED: ${failed[*]}"
  exit 1
fi
