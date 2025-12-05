#!/bin/bash
# Update BADS version number across all files
# Usage: ./scripts/update_version.sh 1.1.4

set -e

if [ -z "$1" ]; then
    echo "Usage: $0 <version>"
    echo "Example: $0 1.1.4"
    exit 1
fi

NEW_VERSION="$1"
DATE=$(date +"%b %-d, %Y")      # e.g., "Dec 5, 2025"
DATE_SHORT=$(date +"%b/%d/%Y")  # e.g., "Dec/05/2025"

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ROOT_DIR="$(dirname "$SCRIPT_DIR")"

echo "Updating to version $NEW_VERSION (date: $DATE)"

# 1. README.md - title
sed -i "s/# Bayesian Adaptive Direct Search (BADS) - v[0-9.]\+/# Bayesian Adaptive Direct Search (BADS) - v$NEW_VERSION/" "$ROOT_DIR/README.md"

# 2. bads.m - function header comment
sed -i "s/%BADS Constrained optimization using Bayesian Adaptive Direct Search (v[0-9.]\+)/%BADS Constrained optimization using Bayesian Adaptive Direct Search (v$NEW_VERSION)/" "$ROOT_DIR/bads.m"

# 3. bads.m - Version in header block
sed -i "s/%   Version: [0-9.]\+/%   Version: $NEW_VERSION/" "$ROOT_DIR/bads.m"

# 4. bads.m - Release date in header block
sed -i "s/%   Release date: .\+/%   Release date: $DATE/" "$ROOT_DIR/bads.m"

# 5. bads.m - bads_version variable
sed -i "s/bads_version = '[0-9.]\+';/bads_version = '$NEW_VERSION';/" "$ROOT_DIR/bads.m"

echo "Done. Updated:"
echo "  - README.md (title)"
echo "  - bads.m (header, version block, variable)"
echo ""
echo "Don't forget to add a changelog entry in bads.m:"
echo "  % $NEW_VERSION ($DATE_SHORT) <description>"
