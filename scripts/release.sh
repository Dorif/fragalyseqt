#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")/.."

usage() {
  cat <<'EOF'
Usage: scripts/release.sh clean|build|check|test-upload|upload

  clean        Remove previous build artifacts.
  build        Build a fresh source distribution and wheel.
  check        Validate metadata and contents with Twine.
  test-upload  Upload dist/* to TestPyPI.
  upload       Upload dist/* to the production PyPI repository.

For uploads, set TWINE_USERNAME=__token__ and TWINE_PASSWORD to the API token.
EOF
}

require_clean_tree() {
  if ! git diff --quiet || ! git diff --cached --quiet; then
    echo "Refusing to publish from a dirty Git tree." >&2
    exit 1
  fi
}

case "${1:-}" in
  clean)
    rm -rf build dist src/*.egg-info
    ;;
  build)
    rm -rf build dist src/*.egg-info
    python scripts/prepare_pypi_readme.py
    python -m build
    ;;
  check)
    python scripts/prepare_pypi_readme.py --check
    python -m twine check --strict dist/*
    ;;
  test-upload)
    require_clean_tree
    python -m twine upload --repository testpypi dist/*
    ;;
  upload)
    require_clean_tree
    python -m twine upload dist/*
    ;;
  *)
    usage >&2
    exit 2
    ;;
esac
