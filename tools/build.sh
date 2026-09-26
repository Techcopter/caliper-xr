#!/usr/bin/env bash
# Validate the extension and build its installable zip into dist/.
# BLENDER overrides the Blender binary (default: the macOS app bundle).
set -euo pipefail
BLENDER="${BLENDER:-/Applications/Blender.app/Contents/MacOS/Blender}"
cd "$(dirname "$0")/.."
mkdir -p dist
"$BLENDER" --command extension validate caliper_xr
"$BLENDER" --command extension build --source-dir caliper_xr --output-dir dist
ls -l dist/*.zip
