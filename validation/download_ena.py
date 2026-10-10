#!/usr/bin/env python3
"""
Download the paired FASTQ files of SRA/ENA runs from the European Nucleotide Archive (ENA) and check every file
against the MD5 sum ENA publishes. Works on Windows, Linux and macOS; needs no SRA Toolkit and no bash.

  python validation/download_ena.py SraRunInfo.csv --out fastq     (a Run Selector CSV with a "Run" column)
  python validation/download_ena.py runs.txt --out fastq           (one run accession per line)
  python validation/download_ena.py SRR13969386 SRR13969385 --out fastq
  python validation/download_ena.py SraRunInfo.csv --out fastq --dry-run     (sizes only, nothing downloaded)
  python validation/download_ena.py runs.tsv --out fastq --max-runs 2        (a pilot: two runs spread over the list; runs.tsv is ENA's run table)

Files are saved as <run>_R1.fastq.gz and <run>_R2.fastq.gz, and the checked sums are written to
<out>/fastq_checksums.md5. Files that are already there with the right checksum are not downloaded again.
Needs access to www.ebi.ac.uk and ftp.sra.ebi.ac.uk.
"""
import argparse
import csv
import hashlib
import http.client
import re
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path
from typing import List

API = ("https://www.ebi.ac.uk/ena/portal/api/filereport?accession={acc}&result=read_run"
       "&fields=run_accession,fastq_ftp,fastq_md5,fastq_bytes&format=tsv")
ACCESSION = re.compile(r"^[SED]RR\d+$")
WAITS = (5, 20, 60)          # seconds to wait before the 2nd, 3rd and 4th attempt to ask ENA about a run
PAUSE = 1                    # multiplies WAITS (the tests set it to 0)


def read_accessions(source: str, ignored: list = None) -> List[str]:
    """Run accessions from a Run Selector CSV (column "Run"), a text file (one per line), or one accession.
    Anything that is not a run accession (a header line such as "acc", say) is skipped and added to `ignored`."""
    path = Path(source)
    if not path.is_file():
        return [source] if ACCESSION.match(source) else []
    text = path.read_text(encoding="utf-8-sig", errors="replace")
    if text.lstrip()[:15].lower().startswith(("<!doctype", "<html")):
        title = re.search(r"<title>(.*?)</title>", text, re.S | re.I)
        raise SystemExit(f"{source} is a web page, not a list of runs" + (f" (its title: {' '.join(title.group(1).split())!r})" if title else "") +
                         ". ENA probably answered with an error when the file was saved. Try again later, and use curl -f so that an "
                         "error is never saved as a file.")
    lines = text.splitlines()
    header = [h.strip().strip('"') for h in lines[0].split(",")] if lines else []
    if "Run" in header:
        accs = [(row.get("Run") or "").strip() for row in csv.DictReader(lines)]
    else:
        accs = [re.split(r"[,\s;]+", ln.strip())[0] for ln in lines if ln.strip() and not ln.lstrip().startswith("#")]
    seen, out = set(), []
    for a in accs:
        if a and not ACCESSION.match(a):
            if ignored is not None:
                ignored.append(a)
            continue
        if a and a not in seen:
            seen.add(a)
            out.append(a)
    return out


def evenly_spaced(items: list, k: int) -> list:
    """k items spread evenly over the list, in order, always including the first and (for k >= 2) the last; all of them when k covers the list."""
    n = len(items)
    if k <= 0 or k >= n:
        return list(items)
    if k == 1:
        return [items[0]]
    return [items[round(i * (n - 1) / (k - 1))] for i in range(k)]


def md5_of(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def ena_failure(accession: str, err: OSError) -> str:
    """What went wrong with a request to ENA, said by kind: a server error is ENA's, a refusal is about the request, anything else is
    the connection."""
    code = getattr(err, "code", None)
    if code is not None and code >= 500:
        return (f"{accession}: ENA answered with a server error (HTTP {code}). The problem is on ENA's side and usually clears within "
                f"minutes: run the same command again later.")
    if code is not None:
        return (f"{accession}: ENA refused the request (HTTP {code}). Check that the accession is correct and that the file you gave "
                f"lists run accessions.")
    return (f"{accession}: could not reach ENA ({err}). Check the internet connection, a proxy or firewall, and that "
            f"www.ebi.ac.uk and ftp.sra.ebi.ac.uk are allowed.")


def paired_files(accession: str, opener=None, attempts: int = len(WAITS) + 1) -> dict:
    """{'1': (url, md5, bytes), '2': (url, md5, bytes)} for a run, as listed by ENA. A server error, no answer or an empty list (ENA
    sometimes answers with one for a run that exists) is tried again after growing pauses; a refusal (a 4xx error) is not."""
    opener = opener or urllib.request.urlopen
    for attempt in range(1, attempts + 1):
        failure, transient = None, True
        try:
            with opener(API.format(acc=accession), timeout=60) as r:
                rows = list(csv.DictReader(r.read().decode().splitlines(), delimiter="\t"))
            if rows and rows[0].get("fastq_ftp"):
                break
            failure = f"{accession}: ENA lists no FASTQ files for this run"
        except OSError as err:                   # includes HTTP and connection errors
            failure = ena_failure(accession, err)
            code = getattr(err, "code", None)
            transient = code is None or code >= 500 or code == 429
        if not transient:
            raise SystemExit(failure)
        if attempt < attempts:
            wait = WAITS[min(attempt - 1, len(WAITS) - 1)] * PAUSE
            print(f"{failure} Trying again in {wait:g} s (attempt {attempt + 1} of {attempts}).")
            time.sleep(wait)
    else:
        raise SystemExit(failure)
    urls, sums = rows[0]["fastq_ftp"].split(";"), rows[0].get("fastq_md5", "").split(";")
    sizes = (rows[0].get("fastq_bytes") or "").split(";")
    pair = {}
    for i, (url, md5) in enumerate(zip(urls, sums)):
        for mate in "12":
            if url.endswith(f"_{mate}.fastq.gz"):
                pair[mate] = (url, md5, int(sizes[i]) if i < len(sizes) and sizes[i].isdigit() else 0)
    if len(pair) != 2:
        raise SystemExit(f"{accession}: ENA does not list a paired FASTQ set (found: {urls})")
    return pair


class ServerAnswer(Exception):
    """The server answered with a short message (an error page) instead of the file. That is not a lost connection, so it is not retried."""


def fetch_resumable(url: str, part: Path, expected_size: int = 0, opener=None) -> None:
    """Download `url` into `part`, continuing from what is already there when the server supports ranges.
    Raises OSError if the transfer ends early, leaving the partial file for the next attempt."""
    opener = opener or urllib.request.urlopen
    start = part.stat().st_size if part.exists() else 0
    if expected_size and start > expected_size:
        part.unlink()
        start = 0
    if not (expected_size and start == expected_size):
        request = urllib.request.Request(url, headers={"Range": f"bytes={start}-"} if start else {})
        try:
            response = opener(request, timeout=60)
        except urllib.error.HTTPError as err:
            if err.code == 416 and start:                  # the server cannot continue from there: the partial file is not a valid beginning (an error page, a changed file)
                part.unlink(missing_ok=True)
                return fetch_resumable(url, part, expected_size, opener)
            raise
        with response:
            resumed = bool(start) and getattr(response, "status", 200) == 206     # 200 means the range was ignored
            with open(part, "ab" if resumed else "wb") as fh:
                while True:
                    data = response.read(1 << 20)
                    if not data:
                        break
                    fh.write(data)
    size = part.stat().st_size
    if expected_size and size < expected_size:
        if size < 4096:
            head = part.read_bytes()[:400]
            if not head.startswith(b"\x1f\x8b"):                       # a gzip file begins with these two bytes; a short answer that does not is a message, not a file
                part.unlink(missing_ok=True)
                raise ServerAnswer(head.decode("utf-8", "replace").strip().replace("\n", " ")[:300] or "an empty answer")
        raise OSError(f"incomplete: got {size} of {expected_size} bytes")


def download(accession: str, out_dir: Path, retries: int = 8, log=print, opener=None, retriever=None,
             file_opener=None) -> List[str]:
    """Download one run's pair into out_dir, verified. Returns 'md5  filename' lines."""
    out_dir.mkdir(parents=True, exist_ok=True)
    records = []
    for mate, (url, md5, size) in paired_files(accession, opener).items():
        target = out_dir / f"{accession}_R{mate}.fastq.gz"
        if target.exists() and md5_of(target) == md5:
            log(f"{accession} R{mate}: already downloaded, checksum ok")
        else:
            part = out_dir / f"{accession}_R{mate}.part"
            interrupted = False
            for attempt in range(1, retries + 1):
                try:
                    if retriever is None:
                        fetch_resumable("https://" + url, part, size, file_opener)
                    else:
                        retriever("https://" + url, part)
                except ServerAnswer as err:
                    part.unlink(missing_ok=True)
                    raise SystemExit(f"{accession} R{mate}: ENA answered with a short message instead of the {size:,} byte file ({err}). This is not a connection problem, so it is not retried. "
                                     f"The run may be missing or restricted at ENA; the same run can be fetched from NCBI (fasterq-dump {accession} --split-files) and put in the raw folder as "
                                     f"{accession}_R1.fastq.gz and {accession}_R2.fastq.gz, or the dataset can be used without it.")
                except (OSError, http.client.HTTPException) as err:
                    interrupted = True
                    log(f"{accession} R{mate}: interrupted ({err}); resuming (attempt {attempt} of {retries})")
                    continue
                interrupted = False
                if md5_of(part) == md5:
                    part.replace(target)
                    log(f"{accession} R{mate}: downloaded, checksum ok")
                    break
                log(f"{accession} R{mate}: checksum mismatch, starting again (attempt {attempt} of {retries})")
                part.unlink(missing_ok=True)
            else:
                if not interrupted:
                    part.unlink(missing_ok=True)
                raise SystemExit(f"{accession} R{mate}: no verified download after {retries} attempts"
                                 + ("; the partial file is kept, so running the same command again resumes it"
                                    if interrupted else "; nothing kept"))
        records.append(f"{md5}  {target.name}")
    return records


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("sources", nargs="+", help="a Run Selector CSV, a text file of accessions, or accessions")
    ap.add_argument("--out", type=Path, required=True, help="folder for the FASTQ files")
    ap.add_argument("--dry-run", action="store_true", help="list the runs and their sizes without downloading")
    ap.add_argument("--max-runs", type=int, default=0, metavar="N", help="use only N runs spread evenly over the list: a pilot, or a small set (default: all)")
    ap.add_argument("--retries", type=int, default=8, help="attempts per file; interrupted transfers resume (default 8)")
    args = ap.parse_args()
    accessions = []
    for s in args.sources:
        skipped = []
        found = read_accessions(s, skipped)
        if skipped:
            print(f"{s}: ignored {', '.join(skipped[:5])} (not run accessions, for example a header line)")
        if not found:
            print(f"warning: no run accessions found in {s}")
        accessions += [a for a in found if a not in accessions]
    if not accessions:
        return 1
    if args.max_runs and len(accessions) > args.max_runs:
        accessions = evenly_spaced(accessions, args.max_runs)
        print(f"--max-runs {args.max_runs}: using {len(accessions)} runs spread evenly over the list: {', '.join(accessions)}")
    if args.dry_run:
        total = 0
        unlisted = []
        for acc in accessions:
            try:
                size = sum(f[2] for f in paired_files(acc).values())
            except SystemExit as err:                       # one run that cannot be listed must not hide the others
                print(err)
                unlisted.append(acc)
                continue
            total += size
            print(f"{acc}: {size / 1e6:,.0f} MB")
        print(f"{len(accessions) - len(unlisted)} of {len(accessions)} run(s) listed, {total / 1e9:,.2f} GB in total. Nothing was downloaded.")
        if unlisted:
            print(f"Could not be listed: {', '.join(unlisted)}.")
        return 1 if unlisted else 0
    records, failed = [], []
    for acc in accessions:
        try:
            records += download(acc, args.out, retries=args.retries)
        except SystemExit as err:                       # one bad run must not stop the others
            failed.append(acc)
            print(err)
    sums = args.out / "fastq_checksums.md5"
    if records:
        sums.write_text("\n".join(records) + "\n", encoding="utf-8")
    done = len(accessions) - len(failed)
    print(f"{done} of {len(accessions)} run(s) downloaded and verified" + (f"; checksums in {sums}" if records else ""))
    if failed:
        print(f"Not completed: {', '.join(failed)}. Run the same command again: verified files are skipped and "
              f"interrupted ones resume.")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
