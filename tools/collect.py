#!/usr/bin/env python3
"""Copy the calculation outputs the explorer needs from the project folder into calculations/.

usage: python3 tools/collect.py --root HOST:/path/to/project     (over ssh, read-only on the remote side)
       python3 tools/collect.py --root /path/to/project          (when run on the cluster itself)

What it writes, for every structure, IRC and path in tools/site_config.py whose source exists:
  calculations/<id>/output.out.gz   the M06-2X output (gzip; host name and working directory redacted)
  calculations/<id>/dlpno.out.gz    the DLPNO-CCSD(T) single point, once it has terminated normally
  calculations/<id>/input.inp       the ORCA input, as echoed at the top of the output
  calculations/<id>/source.txt      where it came from (relative to the project folder) and its state when copied
  calculations/<irc id>/trajectory.xyz, calculations/paths/<id>/trajectory.xyz   path frames with their energies
  calculations/manifest.json        a summary of the above
IRCs are copied only once they have terminated normally; Opt + Freq outputs are copied in any state, so a
structure whose frequencies are still running is included with its electronic energy only.
Nothing on the remote side is written: the files are streamed with tar over ssh.
"""
import argparse
import datetime
import glob as globmod
import gzip
import hashlib
import io
import json
import os
import re
import shlex
import subprocess
import sys
import tarfile
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import site_config as C  # noqa: E402

REPO = os.path.dirname(HERE)
CALC = os.path.join(REPO, "calculations")

# Lines and strings that identify the cluster, the account or the user are redacted before anything is written.
# The account and user names are not written here: set_root() takes them from the --root path at run time.
REDACT_LINES = [(re.compile(r"^(\s*\*\s*Host name:).*$"), r"\1     [redacted]"),
                (re.compile(r"^(\s*\*\s*Working dir\.:).*$"), r"\1  [redacted]")]
REDACT_STRINGS = [(re.compile(r"/nesi/\S*"), "[redacted path]"),
                  (re.compile(r"\S*\.hpc\.nesi\.org\.nz"), "[redacted host]")]
FORBIDDEN_PARTS = [r"/nesi/", r"hpc\.nesi\.org\.nz"]
FORBIDDEN = re.compile("|".join(FORBIDDEN_PARTS))
GENERIC_DIRS = {"nesi", "nobackup", "project", "projects", "home", "scratch", "persistent", "scale_wlg", "Users"}


def set_root(root):
    """Also redact the account and user folders named in the project path (e.g. /<fs>/<area>/<account>/<user>)."""
    global FORBIDDEN
    path = root.split(":", 1)[1] if ":" in root else root
    for part in path.strip("/").split("/"):
        if part and part not in GENERIC_DIRS and len(part) > 3:
            REDACT_STRINGS.append((re.compile(re.escape(part)), "[redacted]"))
            FORBIDDEN_PARTS.append(re.escape(part))
    FORBIDDEN = re.compile("|".join(FORBIDDEN_PARTS))


def redact(text):
    out = []
    for ln in text.split("\n"):
        for rx, rep in REDACT_LINES:
            ln = rx.sub(rep, ln)
        for rx, rep in REDACT_STRINGS:
            ln = rx.sub(rep, ln)
        out.append(ln)
    text = "\n".join(out)
    assert not FORBIDDEN.search(text), "redaction missed something"
    return text


def fetch(root, patterns, dest):
    """Copy every existing file matching the patterns (relative to root) into dest, keeping relative paths."""
    lines = ["set -e", "cd " + shlex.quote(root.split(":", 1)[1] if ":" in root else root),
             "{"]
    for p in patterns:
        # globs are expanded by the remote shell; plain names are quoted
        lit = p if any(ch in p for ch in "*?[") else shlex.quote(p)
        lines.append(f"  for f in {lit}; do [ -f \"$f\" ] && echo \"$f\"; done")
    lines += ["} | tar cf - -T -"]
    script = "\n".join(lines) + "\n"
    if ":" in root:
        host = root.split(":", 1)[0]
        cmd = ["ssh", "-x", "-T", "-o", "BatchMode=yes", "-o", "ConnectTimeout=30", host, "bash -s"]
    else:
        cmd = ["bash", "-s"]
    with tempfile.TemporaryFile() as tf:
        subprocess.run(cmd, input=script.encode(), stdout=tf, check=True)
        tf.seek(0)
        with tarfile.open(fileobj=tf, mode="r:") as tar:
            for m in tar.getmembers():
                if not m.isfile() or m.name.startswith("/") or ".." in m.name.split("/"):
                    continue
                tar.extract(m, dest, filter="data")


def read(path):
    with open(path, errors="replace") as fh:
        return fh.read()


def out_state(text):
    """A short description of an ORCA output's state."""
    normal = "ORCA TERMINATED NORMALLY" in text
    s = []
    if "GEOMETRY OPTIMIZATION" in text or "OptTS" in text or "OPTIMIZATION RUN DONE" in text:
        s.append("optimisation converged" if "THE OPTIMIZATION HAS CONVERGED" in text else "optimisation not converged")
    if "IRC PATH SUMMARY" in text or "THE IRC HAS CONVERGED" in text:
        s.append(f"IRC converged in {text.count('THE IRC HAS CONVERGED')} direction(s)")
    if "VIBRATIONAL FREQUENCIES" in text:
        s.append("frequencies done")
    elif re.search(r"!\s*.*\bFREQ\b", text, re.I) or "%FREQ" in text.upper():
        s.append("frequencies not done yet")
    s.append("terminated normally" if normal else "still running or stopped")
    return ", ".join(s), normal


def input_echo(text):
    lines, on = [], False
    for ln in text.split("\n"):
        if "INPUT FILE" in ln and not on:
            on = True
            continue
        if on:
            if "****END OF INPUT****" in ln:
                break
            m = re.match(r"^\|\s*\d+>\s?(.*)$", ln)
            if m:
                lines.append(m.group(1))
    return "\n".join(lines).rstrip() + "\n"


def write_gz(path, text):
    data = text.encode()
    with open(path, "wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0, compresslevel=9) as gz:
            gz.write(data)
    return hashlib.sha256(data).hexdigest()


def store_output(tmp, src, dest_dir, name, today, note=""):
    p = os.path.join(tmp, src)
    if not os.path.isfile(p):
        return None
    text = redact(read(p))
    state, normal = out_state(text)
    os.makedirs(dest_dir, exist_ok=True)
    sha = write_gz(os.path.join(dest_dir, name), text)
    return dict(file=name, source=src, state=state, normal=normal, sha256=sha, collected=today, note=note)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", required=True, help="HOST:/path/to/project, or a local path on the cluster")
    args = ap.parse_args()
    set_root(args.root)
    today = datetime.date.today().isoformat()

    want = []
    for s in C.SPECIES:
        want += [s["dft"]] + s.get("cc", [])
    for irc in C.IRCS:
        want += [f"{irc['src']}/output.out", f"{irc['src']}/orca_IRC_Full_trj.xyz"]
        for e in irc["ends"]:
            want += [f"{e['src']}/output.out"] + e.get("cc", [])
    for p in C.PATHS:
        want += p["src"] if isinstance(p["src"], list) else [p["src"]]

    manifest = dict(collected=today, root="(project folder)", calculations={})
    with tempfile.TemporaryDirectory() as tmp:
        fetch(args.root, want, tmp)

        def do_calc(cid, dft_src, cc_srcs, extra=None):
            d = os.path.join(CALC, cid)
            rec = dict(id=cid)
            dft = store_output(tmp, dft_src, d, "output.out.gz", today)
            if dft is None:
                print(f"  {cid:<16} not available yet ({dft_src})")
                return None
            rec["dft"] = dft
            with open(os.path.join(d, "input.inp"), "w") as fh:
                fh.write(redact(input_echo(read(os.path.join(tmp, dft_src)))))
            for src in cc_srcs or []:
                p = os.path.join(tmp, src)
                if os.path.isfile(p) and "ORCA TERMINATED NORMALLY" in read(p):
                    rec["cc"] = store_output(tmp, src, d, "dlpno.out.gz", today)
                    break
            if extra:
                rec.update(extra(d))
            lines = [f"Copied on {today} by tools/collect.py from the project folder (paths relative to it).",
                     f"output.out.gz: {dft_src}", f"  state when copied: {dft['state']}"]
            if "cc" in rec:
                lines += [f"dlpno.out.gz: {rec['cc']['source']}", f"  state when copied: {rec['cc']['state']}"]
            lines.append("The host name and working directory in the outputs are redacted.")
            with open(os.path.join(d, "source.txt"), "w") as fh:
                fh.write("\n".join(lines) + "\n")
            manifest["calculations"][cid] = rec
            print(f"  {cid:<16} {dft['state']}" + ("; DLPNO-CCSD(T) done" if "cc" in rec else ""))
            return rec

        print("structures:")
        for s in C.SPECIES:
            do_calc(s["id"], s["dft"], s.get("cc"))
        print("IRCs and their end points:")
        for irc in C.IRCS:
            out = os.path.join(tmp, irc["src"], "output.out")
            if os.path.isfile(out) and "ORCA TERMINATED NORMALLY" in read(out):
                def traj(d, irc=irc):
                    t = os.path.join(tmp, irc["src"], "orca_IRC_Full_trj.xyz")
                    with open(os.path.join(d, "trajectory.xyz"), "w") as fh:
                        fh.write(read(t))
                    return dict(trajectory="trajectory.xyz")
                do_calc(irc["id"], f"{irc['src']}/output.out", None, traj)
            else:
                print(f"  {irc['id']:<16} IRC not finished yet; left out")
            for e in irc["ends"]:
                do_calc(e["id"], f"{e['src']}/output.out", e.get("cc"))
        print("paths:")
        for p in C.PATHS:
            files = p["src"] if isinstance(p["src"], list) else sorted(
                os.path.relpath(f, tmp) for f in globmod.glob(os.path.join(tmp, p["src"])))
            files = [f for f in files if os.path.isfile(os.path.join(tmp, f))]
            if not files:
                print(f"  {p['id']:<16} no frames found")
                continue
            d = os.path.join(CALC, "paths", p["id"])
            os.makedirs(d, exist_ok=True)
            with open(os.path.join(d, "trajectory.xyz"), "w") as fh:
                for f in files:
                    txt = read(os.path.join(tmp, f))
                    fh.write(txt if txt.endswith("\n") else txt + "\n")
            with open(os.path.join(d, "source.txt"), "w") as fh:
                fh.write(f"Copied on {today} by tools/collect.py from the project folder (paths relative to it):\n"
                         + "\n".join(files) + "\n")
            manifest["calculations"]["paths/" + p["id"]] = dict(id=p["id"], frames_from=files, collected=today)
            print(f"  {p['id']:<16} {len(files)} file(s)")

    # keep entries collected earlier whose sources are no longer listed or reachable
    old = os.path.join(CALC, "manifest.json")
    if os.path.isfile(old):
        prev = json.load(open(old)).get("calculations", {})
        for k, v in prev.items():
            manifest["calculations"].setdefault(k, v)
    with open(old, "w") as fh:
        json.dump(manifest, fh, indent=1, sort_keys=True)
    for dirpath, _, files in os.walk(CALC):   # final check: nothing identifying in any text file
        for f in files:
            if f.endswith((".txt", ".inp", ".xyz", ".json")):
                assert not FORBIDDEN.search(read(os.path.join(dirpath, f))), os.path.join(dirpath, f)
    print("written to", os.path.relpath(CALC, os.getcwd()))


if __name__ == "__main__":
    main()
