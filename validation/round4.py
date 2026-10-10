"""Round 4 in one place: every open question of the validation, run in phases with a break between them.

  python validation\\check_build.py --round4 [--quick] [--break-seconds 90] [--phases A,B] [--fresh]

  Phase A  quick analyses of results already on disk (no DADA2; minutes)
  Phase B  simulations with known truth (DADA2; hours)
  Phase C  the standard battery on any dataset registered since the last run (nothing to do when there is none)

Between the phases the run stops for --break-seconds (default 90) and shows what the phase found. During the break:
  Enter or C   continue now          S   skip the next phase          Q   stop here (the report is still written)          P   hold the run until you press Enter
Press nothing and the run CONTINUES by itself when the time is up, so a night's run never waits for you. --break-seconds 0 removes the breaks.
(When the run is not at a keyboard, the breaks simply last their time.)

The steps are those of round4_steps.py, registered in overnight.py as tier 4; their criteria were fixed in TEST_REGISTER.md (section 9) before any run.
It writes ROUND4_REPORT.md in the output folder: per phase, the status of each step, its verdicts, and what each outcome means for the next change.
Run it with `--quick --phases A,B` first (about 25 minutes): the same steps on small slices, to see that every path and tool works (a quick run reads the earlier results of the full run;
phase C is left out because it would only repeat the standard battery on small slices).
"""
import os
import subprocess
import sys
import time
from pathlib import Path
from types import SimpleNamespace

HERE = Path(__file__).resolve().parent

PHASES = [
    ("A", "quick analyses of results already on disk (no DADA2; minutes)", ["split_tiebreak", "missing_variants", "figaro_usable_pairs", "asv_length_check", "figaro_rule_check"]),
    ("B", "simulations with known truth (DADA2; hours)", ["richness_sim", "overmerge_sim"]),
    ("C", "datasets registered since the last run (the standard battery; nothing to do when there is none)", ["battery", "extra_datasets"]),
]

QUESTIONS = {
    "split_tiebreak": "Q28 forward-heavy against balanced split",
    "missing_variants": "Q29 the variant that only one tool recovers",
    "figaro_usable_pairs": "Q30 pairs usable by a tool that needs reads of one length",
    "asv_length_check": "Q23 ASV lengths of the mocks (false joins)",
    "figaro_rule_check": "Q16 FIGARO's maxEE rule reproduced from the reads",
    "richness_sim": "Q0-14 richness against the yield estimate",
    "overmerge_sim": "Q21-27 wrong-length merges and the 12 nt overlap",
    "battery": "the standard battery on the original datasets",
    "extra_datasets": "the standard battery on new datasets",
}

# what a verdict means for the next change (the decisions fixed in TEST_REGISTER.md, section 9)
FOLLOW_UP = {
    ("split_tiebreak", "met"): "supported: add the opt-in `--split-preference forward` to the tool, and test it on the next datasets before changing the default",
    ("split_tiebreak", "NOT MET"): "not supported: keep the balanced tie-break",
    ("figaro_rule_check", "met"): "FIGARO's rule is reproduced: an opt-in `--auto-max-ee fitted` is justified; test it against FIGARO as shipped (merged pairs, and variants on the mocks)",
    ("figaro_rule_check", "NOT MET"): "FIGARO's rule is not reproduced from the reads: no automatic maxEE is built; users can pass FIGARO's own maxEE with --max-ee",
    ("richness_sim", "met"): "the excess grows with richness (or partly): the README says the yield level is an upper bound that rises with the number of taxa",
    ("richness_sim", "NOT MET"): "richness does not explain the excess: the README says the cause is unknown",
    ("overmerge_sim", "met"): "the 12 nt minimum overlap stands in this simulation; say so, and keep it",
    ("overmerge_sim", "NOT MET"): "wrong-length merges or a loss of taxa at 12 nt: consider a longer minimum overlap, and check it on real data before changing the default",
    ("asv_length_check", "met"): "no false-join signature in mock A",
    ("asv_length_check", "NOT MET"): "ASVs of the wrong length in mock A: look at them (the report lists the largest) before anything else",
}


# the verdict of each step that decides something: the text its criterion starts with
DECISION = {
    "split_tiebreak": "a forward-heavy split",
    "figaro_rule_check": "the rule 'mean expected errors",
    "richness_sim": "the excess of the yield estimate grows",
    "overmerge_sim": "the 12 nt minimum overlap stands",
    "asv_length_check": "mock A, AmpliTrim's default setting",
}


def default_key_reader():
    """A function returning the key pressed since the last call ('' for none; Enter is '\\n'). Windows: msvcrt; elsewhere a terminal; otherwise no keys."""
    if os.name == "nt":
        try:
            import msvcrt

            def read():
                if msvcrt.kbhit():
                    ch = msvcrt.getwch()
                    return "\n" if ch in ("\r", "\n") else ch.lower()
                return ""
            return read
        except ImportError:
            return lambda: ""
    if sys.stdin is not None and sys.stdin.isatty():
        import select

        def read():
            ready, _, _ = select.select([sys.stdin], [], [], 0)
            if ready:
                line = sys.stdin.readline()
                return line.strip()[:1].lower() or "\n"
            return ""
        return read
    return lambda: ""


def pause_between(seconds: float, read_key, now=time.monotonic, sleep=time.sleep, say=print) -> str:
    """Wait for a key up to `seconds`: 'continue' (Enter, C, or no key at all), 'skip' (S) or 'quit' (Q). P holds the run until Enter, S or Q."""
    if seconds <= 0:
        return "continue"
    deadline, last = now() + seconds, None
    while True:
        remaining = deadline - now()
        if remaining <= 0:
            return "continue"
        if int(remaining) // 10 != last:
            last = int(remaining) // 10
            say(f"   continuing in {int(remaining)} s   (Enter = now, S = skip the next phase, Q = stop here, P = pause)")
        key = read_key()
        if key in ("\n", "c"):
            return "continue"
        if key == "s":
            return "skip"
        if key == "q":
            return "quit"
        if key == "p":
            say("   paused: Enter = continue, S = skip the next phase, Q = stop here")
            while True:
                key = read_key()
                if key in ("\n", "c"):
                    return "continue"
                if key == "s":
                    return "skip"
                if key == "q":
                    return "quit"
                sleep(0.2)
        sleep(0.2)


def read_result(out: Path, step: str):
    import json
    f = Path(out) / step / "result.json"
    if not f.is_file():
        return None
    try:
        return json.loads(f.read_text(encoding="utf-8"))
    except ValueError:
        return None


def step_lines(out: Path, step: str) -> list:
    """What one step found, for the screen and for the report."""
    res = read_result(out, step)
    label = QUESTIONS.get(step, step)
    if res is None:
        return [f"- {step} ({label}): no result"]
    head = f"- {step} ({label}): {res.get('status', '?')}" + (f" in {res.get('seconds', 0) / 60:.1f} min" if res.get("status") == "ok" else "")
    if res.get("status") != "ok":
        return [head + (f": {res['error']}" if res.get("error") else "")]
    lines = [head]
    for v in res.get("verdicts", []):
        lines.append(f"    {v.get('outcome', '?'):8s} {v.get('criterion', '')}: {v.get('detail', '')}")
        prefix = DECISION.get(step)
        action = FOLLOW_UP.get((step, v.get("outcome")))
        if action and prefix and v.get("criterion", "").startswith(prefix):          # only the decision verdict of a step carries a follow-up
            lines.append(f"             -> {action}")
    return lines


def write_report(out: Path, phases_done: list, started: str, notes: list) -> Path:
    L = ["# AmpliTrim round 4: report", "", f"Started {started}; written {time.strftime('%Y-%m-%d %H:%M')}.", ""]
    if notes:
        L += ["## What happened between the phases", ""] + [f"- {n}" for n in notes] + [""]
    for key, title, code, lines in phases_done:
        L += [f"## Phase {key}: {title}", "", f"overnight.py exit code {code} (0 = no step failed).", ""] + lines + [""]
    L += ["The criteria are in TEST_REGISTER.md, section 9; each step's full tables are in its folder next to this file (OVERNIGHT_SUMMARY.md has them all).", ""]
    path = Path(out) / "ROUND4_REPORT.md"
    Path(out).mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(L), encoding="utf-8")
    return path


def run_round4(cfg: dict, args, runner=None, read_key=None, sleep=time.sleep, now=time.monotonic, say=print) -> int:
    """The phases in order, with a break after each; the run goes on when nothing is pressed."""
    import check_build
    runner = runner or subprocess.call
    out = Path(os.path.expandvars(cfg["out"]) + ("_quick" if getattr(args, "quick", False) else ""))
    wanted = {p.strip().upper() for p in (getattr(args, "phases", None) or "A,B,C").split(",") if p.strip()}
    phases = [p for p in PHASES if p[0] in wanted]
    if not phases:
        raise SystemExit("--phases: choose from A, B, C (for example A,B)")
    break_s = getattr(args, "break_seconds", 90)
    read_key = read_key or default_key_reader()
    started, done, notes, skipped, worst = time.strftime("%Y-%m-%d %H:%M"), [], [], set(), 0
    with check_build.keep_awake():
        for i, (key, title, steps) in enumerate(phases):
            if key in skipped:
                notes.append(f"phase {key} skipped at your request")
                continue
            say(f"\n==================== Phase {key}: {title}\n")
            ns = SimpleNamespace(only=",".join(steps), skip=None, quick=getattr(args, "quick", False), fresh=getattr(args, "fresh", False))
            cmd = [sys.executable, str(HERE / "overnight.py")] + check_build.overnight_args(cfg, ns)
            say("$ " + " ".join(cmd) + "\n")
            code = runner(cmd)
            worst = max(worst, code or 0)
            lines = []
            for s in steps:
                lines += step_lines(out, s)
            say("\n--- what phase " + key + " found\n" + "\n".join(lines))
            done.append((key, title, code, lines))
            if i < len(phases) - 1 and break_s > 0:
                nxt = phases[i + 1]
                say(f"\n--- break before phase {nxt[0]} ({nxt[1]})")
                action = pause_between(break_s, read_key, now=now, sleep=sleep, say=say)
                if action == "quit":
                    notes.append(f"stopped at your request after phase {key}")
                    say("stopped at your request")
                    break
                if action == "skip":
                    skipped.add(nxt[0])
                    say(f"phase {nxt[0]} will be skipped")
    path = write_report(out, done, started, notes)
    say(f"\nreport: {path}")
    return worst
