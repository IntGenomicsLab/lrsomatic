#!/usr/bin/env python3
"""Write <sample>.verdict_status.tsv: whether ClairS-TO's Verdict tagged the calls, and from which purity.

One header line and one row: status (applied, skipped, not_run), purity, ploidy, purity_source, message.

--status sets the status directly. Otherwise it comes from Verdict's own tagging logs (--verdict_log): none of them
present means Verdict never ran (not_run); a line saying "not applying verdict tagging" means it ran and tagged nothing
(skipped); anything else means it tagged (applied). The message holds the log lines that explain a skip or a disabled
Verdict, from --verdict_log and --run_log.
"""
import argparse
import os
import re

REASON = re.compile(
    r"not applying verdict tagging|Could not find an optimal purity|no purity/ploidy solution|"
    r"Disabling Verdict|disable_verdict option",
    re.IGNORECASE,
)
SKIPPED = re.compile(r"not applying verdict tagging", re.IGNORECASE)


def read_lines(paths):
    lines = []
    for path in paths:
        if path and os.path.isfile(path):
            with open(path, errors="replace") as fp:
                lines.extend(line.strip() for line in fp)
    return lines


def purity_ploidy(path):
    """Purity and ploidy from a Verdict-layout table (Sample, Purity, Ploidy, ...); NA when absent."""
    if not path or not os.path.isfile(path):
        return "NA", "NA"
    with open(path) as fp:
        rows = [line.rstrip("\n").split("\t") for line in fp if line.strip()]
    if len(rows) < 2:
        return "NA", "NA"
    values = dict(zip(rows[0], rows[1]))
    return values.get("Purity", "NA") or "NA", values.get("Ploidy", "NA") or "NA"


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--out", required=True)
    parser.add_argument("--purity_source", required=True, help="e.g. ascat, clairsto_internal")
    parser.add_argument("--purity_ploidy", help="the <sample>_Tumor_Purity_Ploidy.txt Verdict read, if any")
    parser.add_argument("--verdict_log", nargs="*", default=[], help="logs of Verdict's tagging step; absent if it did not run")
    parser.add_argument("--run_log", nargs="*", default=[], help="other logs searched for why Verdict was disabled")
    parser.add_argument("--status", choices=["applied", "skipped", "not_run"])
    args = parser.parse_args()

    verdict_lines = read_lines(args.verdict_log)
    if args.status:
        status = args.status
    elif not any(os.path.isfile(p) for p in args.verdict_log):
        status = "not_run"
    elif any(SKIPPED.search(line) for line in verdict_lines):
        status = "skipped"
    else:
        status = "applied"

    message = []
    if status != "applied":
        for line in verdict_lines + read_lines(args.run_log):
            line = re.sub(r"^\[[A-Z]+\]\s*", "", " ".join(line.split()))
            if REASON.search(line) and line not in message:
                message.append(line)

    purity, ploidy = purity_ploidy(args.purity_ploidy)
    source = "none" if status == "not_run" else args.purity_source
    with open(args.out, "w") as out:
        out.write("status\tpurity\tploidy\tpurity_source\tmessage\n")
        out.write("\t".join([status, purity, ploidy, source, " ".join(message)]) + "\n")
    print(f"[INFO] Verdict {status}: purity {purity}, ploidy {ploidy} ({source})")


if __name__ == "__main__":
    main()
