#!/usr/bin/env python3
"""
Install the Windows build of vsearch into the active Python environment (for conda: <env>\\Scripts, which is on
PATH whenever the environment is active). vsearch.exe needs the two DLLs that come with it; they are installed
beside it.

Integrity: a checksum you give with --sha256 is used if present, otherwise the SHA256SUMS published with the
release. If a published checksum exists and does not match the download, nothing is installed unless you pass
--trust-download. If the release publishes no checksums (all but v2.32.0 at the time of writing), the download
is installed with a warning. The SHA256 of what was installed is always printed.

    python validation/install_vsearch_windows.py [--version 2.32.0] [--dest DIR]
"""
import argparse
import hashlib
import io
import platform
import sys
import urllib.error
import urllib.request
import zipfile
from pathlib import Path

BASE = "https://github.com/torognes/vsearch/releases/download/v{v}/{name}"


def fetch(url: str) -> bytes:
    with urllib.request.urlopen(url, timeout=120) as r:
        return r.read()


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--version", default="2.32.0")
    ap.add_argument("--dest", type=Path, default=Path(sys.prefix) / "Scripts")
    ap.add_argument("--sha256", help="the expected SHA256 of the zip, from a source you trust")
    ap.add_argument("--trust-download", action="store_true", help="install even if the published checksum does not match")
    ap.add_argument("--force", action="store_true", help="run on a non-Windows system (for testing)")
    args = ap.parse_args()
    if platform.system() != "Windows" and not args.force:
        print("This installs the Windows build. On Linux or macOS use: conda install -c bioconda vsearch")
        return 1
    name = f"vsearch-{args.version}-win-x86_64.zip"
    print(f"downloading {name} ...")
    data = fetch(BASE.format(v=args.version, name=name))
    actual = hashlib.sha256(data).hexdigest()
    expected, source = args.sha256, "supplied"
    if expected is None:
        try:
            sums = fetch(BASE.format(v=args.version, name="SHA256SUMS")).decode()
            expected = next((ln.split()[0] for ln in sums.splitlines() if ln.strip().endswith(name)), None)
            source = "published"
        except urllib.error.HTTPError:
            expected = None                                   # this release publishes no checksums
    if expected is None:
        print(f"warning: no checksum is published for {name}, so the download cannot be verified.")
    elif actual.lower() == expected.lower():
        print(f"checksum verified ({source})")
    elif args.trust_download:
        print(f"warning: the {source} checksum ({expected}) does not match the download; installing because "
              f"--trust-download was given.")
    else:
        print(f"the {source} checksum does not match the download. Nothing was installed.\n"
              f"  expected {expected}\n  actual   {actual}\n"
              f"Use --sha256 with a hash you trust, --trust-download to accept this file, or --version 2.31.0.")
        return 1
    print(f"sha256 of the download: {actual}")
    args.dest.mkdir(parents=True, exist_ok=True)
    installed = []
    with zipfile.ZipFile(io.BytesIO(data)) as z:
        for info in z.infolist():
            parts = Path(info.filename).parts
            if len(parts) >= 3 and parts[-2] == "bin" and not info.is_dir():
                (args.dest / parts[-1]).write_bytes(z.read(info))
                installed.append(parts[-1])
    if "vsearch.exe" not in installed:
        print("vsearch.exe was not found in the download.")
        return 1
    print(f"installed {', '.join(sorted(installed))} to {args.dest}\nCheck it with: vsearch --version")
    return 0


if __name__ == "__main__":
    sys.exit(main())
