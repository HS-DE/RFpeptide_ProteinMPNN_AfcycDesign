#!/bin/bash
set -euo pipefail

SOURCE_DIR="${COLABDESIGN_GAMMA_SOURCE:-$HOME/fga_model_envs/sources/ColabDesign-gamma-stage5}"
PINNED_COMMIT="5ab4efaba2321a6c3c314b82d2fff8e0241f5c2d"
COMMIT_MARKER="$SOURCE_DIR/.stage5_colabdesign_commit"

if [[ -f "$COMMIT_MARKER" ]]; then
  current_commit="$(tr -d '[:space:]' < "$COMMIT_MARKER")"
  if [[ "$current_commit" != "$PINNED_COMMIT" ]]; then
    echo "FAIL: existing gamma source is at $current_commit, expected $PINNED_COMMIT" >&2
    exit 2
  fi
  if [[ -f "$SOURCE_DIR/colabdesign/af/contrib/predict.py" ]] && [[ -f "$SOURCE_DIR/colabdesign/af/contrib/cyclic.py" ]]; then
    echo "PASS: pinned ColabDesign gamma source already exists."
    exit 0
  fi
  echo "FAIL: commit marker exists but required gamma prediction files are missing." >&2
  exit 2
fi

if [[ -e "$SOURCE_DIR" ]]; then
  echo "FAIL: source path exists without a valid Stage 5 commit marker: $SOURCE_DIR" >&2
  exit 2
fi

mkdir -p "$(dirname "$SOURCE_DIR")"
TMP_DIR="$(mktemp -d)"
trap 'rm -rf "$TMP_DIR"' EXIT
curl -fL --retry 5 --retry-delay 2 \
  "https://codeload.github.com/sokrypton/ColabDesign/tar.gz/$PINNED_COMMIT" \
  -o "$TMP_DIR/colabdesign.tar.gz"
tar -xzf "$TMP_DIR/colabdesign.tar.gz" -C "$TMP_DIR"
EXTRACTED_DIR="$(find "$TMP_DIR" -mindepth 1 -maxdepth 1 -type d -name 'ColabDesign-*' | head -n 1)"
if [[ -z "$EXTRACTED_DIR" ]]; then
  echo "FAIL: downloaded archive did not contain the expected ColabDesign source directory." >&2
  exit 2
fi
mv "$EXTRACTED_DIR" "$SOURCE_DIR"
printf '%s
' "$PINNED_COMMIT" > "$COMMIT_MARKER"
echo "Prepared pinned ColabDesign gamma source at $SOURCE_DIR"
