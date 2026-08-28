#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "$SCRIPT_DIR/.." && pwd)"
ENGINE="${CONTAINER_ENGINE:-docker}"
OUT_DIR="${OUT_DIR:-$REPO_ROOT/dist/packages}"
TARGETS=("${@:-all}")

VERSION="$({ cd "$REPO_ROOT"; python3 - <<'PY'
import re
from pathlib import Path
text = Path('src/fragalyseqt/__init__.py').read_text(encoding='utf-8')
print(re.search(r'^__version__ = "([^"]+)"$', text, re.MULTILINE).group(1))
PY
})"

mkdir -p "$OUT_DIR"/{ubuntu,rpm,altlinux}

container() {
  local image=$1 script=$2
  "$ENGINE" run --rm \
    -e VERSION="$VERSION" \
    -e SOURCE_DATE_EPOCH="${SOURCE_DATE_EPOCH:-1778112000}" \
    -v "$REPO_ROOT:/src:ro" \
    -v "$OUT_DIR:/out" \
    "$image" bash -euxo pipefail -c "$script"
}

build_ubuntu() {
  container ubuntu:22.04 '
    export DEBIAN_FRONTEND=noninteractive
    apt-get update
    apt-get install -y --no-install-recommends build-essential debhelper devscripts fakeroot python3
    mkdir -p /work
    cp -r /src/src /src/debian /src/packaging /work/
    cp /src/README.md /src/COPYING /src/LICENSE /work/
    cd /work
    rm -rf debian/fragalyseqt debian/.debhelper debian/debhelper-build-stamp debian/files debian/*.substvars
    dpkg-buildpackage -us -uc -b
    cp /fragalyseqt_*_all.deb /out/ubuntu/
  '
}

rpm_source_tree() {
  cat <<'SCRIPT'
    mkdir -p "/work/fragalyseqt-$VERSION" /root/rpmbuild/{BUILD,RPMS,SOURCES,SPECS,SRPMS}
    cp -r /src/src /src/packaging "/work/fragalyseqt-$VERSION/"
    cp /src/README.md /src/COPYING /src/LICENSE "/work/fragalyseqt-$VERSION/"
    tar --sort=name --mtime="@$SOURCE_DATE_EPOCH" --owner=0 --group=0 --numeric-owner \
      -C /work -czf "/root/rpmbuild/SOURCES/fragalyseqt-$VERSION.tar.gz" \
      "fragalyseqt-$VERSION"
    cp "/work/fragalyseqt-$VERSION/packaging/fragalyseqt.spec" /root/rpmbuild/SPECS/
SCRIPT
}

# The universal RPM: built on Fedora, installs on every RPM distribution
# except ALT (Fedora, RHEL and derivatives, openSUSE Tumbleweed/Slowroll).
# No dist suffix — the package is not Fedora-specific.
build_rpm() {
  container fedora:latest "
    dnf install -y rpm-build tar gzip python3
$(rpm_source_tree)
    rpmbuild -bb --define '_topdir /root/rpmbuild' --define 'dist %{nil}' /root/rpmbuild/SPECS/fragalyseqt.spec
    cp /root/rpmbuild/RPMS/noarch/*.rpm /out/rpm/
  "
}

# The ALT RPM: ALT uses python3-module-* naming, has no Biopython package and
# a dependency generator that cannot process this package, so it needs a build
# of its own. Installs on ALT only.
build_altlinux() {
  container alt:p11 "
    apt-get update
    apt-get install -y rpm-build tar gzip python3 python3-module-pip shadow-utils
$(rpm_source_tree)
    useradd -m builder
    mv /root/rpmbuild /home/builder/rpmbuild
    chown -R builder:builder /home/builder/rpmbuild
    runuser -u builder -- rpmbuild -bb --define '_topdir /home/builder/rpmbuild' --define 'dist .altlinux' /home/builder/rpmbuild/SPECS/fragalyseqt.spec
    cp /home/builder/rpmbuild/RPMS/noarch/*.rpm /out/altlinux/
  "
}

run_target() {
  case "$1" in
    all) build_ubuntu; build_rpm; build_altlinux ;;
    ubuntu|deb) build_ubuntu ;;
    rpm|fedora|opensuse|suse) build_rpm ;;
    altlinux|alt) build_altlinux ;;
    *) echo "Unknown target: $1" >&2; exit 2 ;;
  esac
}

for target in "${TARGETS[@]}"; do
  run_target "$target"
done

# -print rather than -printf: the latter is a GNU extension that BSD find
# (macOS) does not have, and under `set -e` it aborted the script *after* a
# successful build, making a finished set of packages look like a failure.
find "$OUT_DIR" -type f \( -name '*.deb' -o -name '*.rpm' \) -print | sort
