#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "$SCRIPT_DIR/.." && pwd)"
ENGINE="${CONTAINER_ENGINE:-docker}"
PKG_DIR="${PKG_DIR:-$REPO_ROOT/dist/packages}"
TARGETS=("${@:-all}")

container() {
  local image=$1 package_dir=$2 script=$3
  "$ENGINE" run --rm \
    -v "$PKG_DIR/$package_dir:/packages:ro" \
    "$image" bash -euxo pipefail -c "$script"
}

smoke_gui() {
  cat <<'SCRIPT'
    set +e
    QT_QPA_PLATFORM=offscreen timeout 8s fragalyseqt >/tmp/fragalyseqt.log 2>&1
    rc=$?
    set -e
    if [ "$rc" -ne 124 ]; then
      cat /tmp/fragalyseqt.log
      exit "$rc"
    fi
    ! grep -q Traceback /tmp/fragalyseqt.log
SCRIPT
}

# The DEB is built on 22.04 and must install on 22.04 and every later release.
test_ubuntu() {
  local image=${1:-ubuntu:22.04}
  container "$image" ubuntu "
    export DEBIAN_FRONTEND=noninteractive
    apt-get update
    apt-get install -y --no-install-recommends /packages/*.deb
$(smoke_gui)
    dpkg-query -W fragalyseqt
  "
}

# The SAME universal RPM is verified on every non-ALT RPM distribution.
test_fedora() {
  container fedora:latest rpm "
    dnf install -y /packages/*.rpm
$(smoke_gui)
    rpm -q fragalyseqt
  "
}

test_opensuse() {
  container opensuse/tumbleweed:latest rpm "
    zypper -n --gpg-auto-import-keys refresh
    zypper -n --no-gpg-checks install /packages/*.rpm
$(smoke_gui)
    rpm -q fragalyseqt
  "
}

# Slowroll is not published as a container image, so it is reproduced by
# pointing a Tumbleweed base at the Slowroll repositories.
test_slowroll() {
  container opensuse/tumbleweed:latest rpm "
    zypper -n rr repo-oss repo-non-oss repo-update >/dev/null 2>&1 || true
    zypper -n --gpg-auto-import-keys ar -f https://download.opensuse.org/slowroll/repo/oss/ slowroll-oss
    zypper -n --gpg-auto-import-keys --no-gpg-checks refresh
    zypper -n --no-gpg-checks install /packages/*.rpm
$(smoke_gui)
    rpm -q fragalyseqt
  "
}

test_altlinux() {
  container alt:p11 altlinux "
    apt-get update
    apt-get install -y /packages/*.rpm
$(smoke_gui)
    rpm -q fragalyseqt
    echo '--- Biopython must be present after %post ---'
    python3 -c 'import Bio; print(\"Bio\", Bio.__version__)'
  "
}

for target in "${TARGETS[@]}"; do
  case "$target" in
    all)
      test_ubuntu ubuntu:22.04
      test_ubuntu ubuntu:24.04
      test_fedora
      test_opensuse
      test_slowroll
      test_altlinux
      ;;
    ubuntu|deb) test_ubuntu ubuntu:22.04; test_ubuntu ubuntu:24.04 ;;
    ubuntu22.04) test_ubuntu ubuntu:22.04 ;;
    ubuntu24.04) test_ubuntu ubuntu:24.04 ;;
    rpm|fedora) test_fedora ;;
    opensuse|suse|tumbleweed) test_opensuse ;;
    slowroll) test_slowroll ;;
    alt|altlinux) test_altlinux ;;
    *) echo "Unknown target: $target" >&2; exit 2 ;;
  esac
done
