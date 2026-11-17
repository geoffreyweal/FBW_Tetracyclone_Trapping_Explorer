#!/usr/bin/env python3
"""Build the explorer's data, thermochemistry and tables from calculations/ (see tools/collect.py).

usage: python3 tools/build_site.py          (needs GoodVibes 4.3.0: pip install "goodvibes==4.3.0")

Method, as in Part II of the Supplementary Information:
- Structures: M06-2X/def2-TZVP D3(0), RIJCOSX (def2/J), DefGrid2, TightSCF/TightOPT, CPCM(CH2Cl2: eps 9.08, n 1.424),
  ORCA 6.1.0. Every input is checked against this before anything else is done.
- Thermochemistry: GoodVibes 4.3.0 at 313.15 K from the M06-2X harmonic frequencies, Grimme quasi-RRHO entropy with a
  100 cm-1 cut-off, no quasi-harmonic enthalpy correction, no frequency scaling (-v 1.0). G is GoodVibes' qh-G.
  1 mol/L standard state (-c 1.0, i.e. RT ln(Vm/V°) = 8.45 kJ/mol added to each species); a second run at 1 atm gives
  the 1 atm values.
- DLPNO-CCSD(T)/def2-TZVPP single points, where finished: composite G = Eel(DLPNO-CCSD(T)) + [G - Eel](M06-2X), from
  GoodVibes --spc. Relative energies stay at the M06-2X level until every structure has its DLPNO-CCSD(T) energy.
- Relative energies in kJ/mol with 1 Eh = 2625.50 kJ/mol, relative to the separated reactants of each step.

Writes data/data.js, data/paths.js, xyz/, thermochemistry/ and README.md. index.html and assets/ are hand-written.
"""
import csv
import gzip
import importlib.metadata
import io
import json
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
from datetime import date

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import site_config as C  # noqa: E402

REPO = os.path.dirname(HERE)
CALC = os.path.join(REPO, "calculations")
THERMO = os.path.join(REPO, "thermochemistry")
KJ = 2625.50                       # 1 Eh in kJ/mol, as in the SI
TEMP = 313.15
GV_VERSION = "4.3.0"
GV_ARGS = ["--temp", "313.15", "-v", "1.0", "--imag", "--check"]
SITE_URL = "https://geoffreyweal.github.io/FBW_Tetracyclone_Trapping_Explorer/"


# ------------------------------------------------------------------------------------------------ reading outputs
def read_gz(path):
    with gzip.open(path, "rt", errors="replace") as fh:
        return fh.read()


def last_energy(text):
    e = None
    for m in re.finditer(r"^FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", text, re.M):
        e = float(m.group(1))
    return e


def last_geometry(text):
    blocks = list(re.finditer(r"CARTESIAN COORDINATES \(ANGSTROEM\)\n-+\n((?:\s*[A-Z][a-z]?\s+-?\d+\.\d+\s+-?\d+\.\d+\s+-?\d+\.\d+\n)+)", text))
    if not blocks:
        return None
    return [(t[0], float(t[1]), float(t[2]), float(t[3])) for t in (ln.split() for ln in blocks[-1].group(1).strip().split("\n"))]


def last_freqs(text):
    i = text.rfind("VIBRATIONAL FREQUENCIES")
    if i < 0:
        return None
    out = []
    for m in re.finditer(r"^\s*\d+:\s+(-?\d+\.\d+) cm\*\*-1", text[i:i + 200000], re.M):
        out.append(float(m.group(1)))
        if len(out) > 3 * 400:
            break
    return out


def input_lines(text):
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
    return lines


def method_check(text, kind):
    """Compare an output's input and header with the SI's method. Returns (ok, details)."""
    inp = "\n".join(re.sub(r"#.*$", "", ln) for ln in input_lines(text))
    up = inp.upper()
    ver = re.search(r"Program Version (\S+)", text)
    d = dict(version=ver.group(1) if ver else "?")
    checks = []
    if kind == "dlpno":
        checks = [("DLPNO-CCSD(T)/def2-TZVPP, def2-TZVPP/C, TightPNO",
                   all(k in up for k in ("DLPNO-CCSD(T)", "DEF2-TZVPP", "DEF2-TZVPP/C", "TIGHTPNO"))),
                  ("RIJCOSX def2/J", "RIJCOSX" in up and "DEF2/J" in up),
                  ("TightSCF", "TIGHTSCF" in up), ("DefGrid2", "DEFGRID2" in up)]
    else:
        checks = [("M06-2X/def2-TZVP", "M062X" in up and re.search(r"\bDEF2-TZVP\b", up) is not None),
                  ("D3(0)", "D3ZERO" in up),
                  # RIJCOSX is ORCA's default for hybrids: check the output really built COSX grids and used def2/J
                  ("RIJCOSX with def2/J", "COSX GRID GENERATION" in text and "NORIJCOSX" not in up
                   and "Your calculation utilizes the auxiliary basis: def2/J" in text),
                  ("TightSCF", "TIGHTSCF" in up), ("DefGrid2", "DEFGRID2" in up)]
        if kind == "opt":
            checks.append(("TightOPT", "TIGHTOPT" in up))
    checks.append(("CPCM eps 9.08, n 1.424", re.search(r"%CPCM\s+EPSILON\s+9\.08\s+REFRAC\s+1\.424", up) is not None))
    checks.append(("ORCA 6.1.0", d["version"].startswith("6.1.0")))
    d["method_line"] = " ".join(ln.strip() for ln in input_lines(text) if ln.strip().startswith("!"))
    d["freq"] = ("numerical" if "NUMFREQ" in up else "analytic" if re.search(r"\bFREQ\b", up) else
                 "IRC" if re.search(r"\bIRC\b", up) else "none")
    d["failed"] = [name for name, ok in checks if not ok]
    return not d["failed"], d


def dlpno_info(text):
    t1 = re.search(r"T1 diagnostic\s+\.\.\.\s+(\d+\.\d+)", text)
    return dict(E=last_energy(text), T1=float(t1.group(1)) if t1 else None,
                normal="ORCA TERMINATED NORMALLY" in text)


def read_frames(path):
    lines = open(path).read().split("\n")
    out, i = [], 0
    while i < len(lines) and lines[i].strip():
        n = int(lines[i].split()[0])
        atoms = [(t[0], float(t[1]), float(t[2]), float(t[3])) for t in (ln.split() for ln in lines[i + 2:i + 2 + n])]
        out.append((lines[i + 1].strip(), atoms))
        i += n + 2
    return out


def comment_energy(comment):
    m = re.search(r"E\s+(-?\d+\.\d+)", comment)
    return float(m.group(1)) if m else None


# ------------------------------------------------------------------------------------------------ geometry
def dist(atoms, i, j):
    a, b = atoms[i], atoms[j]
    return math.sqrt((a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2 + (a[3] - b[3]) ** 2)


def kabsch_rmsd(A, B):
    import numpy as np
    P = np.array([a[1:] for a in A]); Q = np.array([b[1:] for b in B])
    P -= P.mean(0); Q -= Q.mean(0)
    U, S, Vt = np.linalg.svd(P.T @ Q)
    dd = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, dd]) @ U.T
    return float(np.sqrt(((P @ R.T - Q) ** 2).sum(1).mean()))


def distance_spectrum_rms(A, B):
    """Permutation-invariant comparison for structures whose atom order differs (sorted distances per element pair)."""
    import numpy as np
    if sorted(a[0] for a in A) != sorted(b[0] for b in B):
        return None
    def spec(X):
        out = {}
        for i in range(len(X)):
            for j in range(i + 1, len(X)):
                out.setdefault(tuple(sorted((X[i][0], X[j][0]))), []).append(dist(X, i, j))
        return {k: np.sort(v) for k, v in out.items()}
    a, b = spec(A), spec(B)
    diffs = np.concatenate([a[k] - b[k] for k in a])
    return float(np.sqrt((diffs ** 2).mean()))


def compare(A, B):
    """(measure, value, heavy-atom value) for two structures."""
    if [a[0] for a in A] == [b[0] for b in B]:
        heavy = [i for i, a in enumerate(A) if a[0] != "H"]
        return "RMSD", kabsch_rmsd(A, B), kabsch_rmsd([A[i] for i in heavy], [B[i] for i in heavy])
    return "distance RMS (atom order differs)", distance_spectrum_rms(A, B), None


def xyz_text(atoms, comment):
    return f"{len(atoms)}\n{comment}\n" + "\n".join(f"{a[0]:<2} {a[1]:14.8f} {a[2]:14.8f} {a[3]:14.8f}" for a in atoms) + "\n"


# ------------------------------------------------------------------------------------------------ GoodVibes
def check_goodvibes():
    try:
        v = importlib.metadata.version("goodvibes")
    except importlib.metadata.PackageNotFoundError:
        sys.exit(f"GoodVibes is not installed: pip install 'goodvibes=={GV_VERSION}'")
    if v != GV_VERSION:
        sys.exit(f"GoodVibes {v} found; the SI uses {GV_VERSION}: pip install 'goodvibes=={GV_VERSION}'")
    return v


def run_goodvibes(ids, tag, conc, spc=False):
    """Run GoodVibes on calculations/<id>/output.out.gz (and dlpno.out.gz with --spc). Returns {id: result}."""
    with tempfile.TemporaryDirectory() as tmp:
        for cid in ids:
            with open(os.path.join(tmp, f"{cid}.out"), "w") as fh:
                fh.write(read_gz(os.path.join(CALC, cid, "output.out.gz")))
            if spc:
                with open(os.path.join(tmp, f"{cid}_CC.out"), "w") as fh:
                    fh.write(read_gz(os.path.join(CALC, cid, "dlpno.out.gz")))
        args = GV_ARGS + (["-c", "1.0"] if conc else []) + (["--spc", "CC"] if spc else [])
        if spc:   # GoodVibes 4.3.0's --check crashes with --spc (TypeError in validation.check_files); the DFT runs check
            args = [a for a in args if a != "--check"]
        cmd = [sys.executable, "-I", "-m", "goodvibes", *args, "--output", tag, "--json", f"{tag}.json",
               *[f"{cid}.out" for cid in ids]]
        r = subprocess.run(cmd, cwd=tmp, env=dict(os.environ, COLUMNS="250"), capture_output=True, text=True)
        if r.returncode != 0 or not os.path.isfile(os.path.join(tmp, f"{tag}.json")):
            sys.exit(f"GoodVibes failed ({tag}):\n{r.stdout}\n{r.stderr}")
        js = json.load(open(os.path.join(tmp, f"{tag}.json")))
        dat = open(os.path.join(tmp, f"GoodVibes_{tag}.dat")).read()
    for res in js["results"]:   # keep file names only
        for d in (res, res.get("qcdata", {})):
            for k in ("file", "sp_file"):
                if d.get(k):
                    d[k] = os.path.basename(d[k])
    js.pop("generated_at", None)
    os.makedirs(THERMO, exist_ok=True)
    with open(os.path.join(THERMO, f"GoodVibes_{tag}.dat"), "w") as fh:
        fh.write(f"# goodvibes {' '.join(args)} <structure>.out ...   (run by tools/build_site.py)\n" + dat)
    with open(os.path.join(THERMO, f"goodvibes_{tag}.json"), "w") as fh:
        json.dump(js, fh, indent=1)
    out = {}
    for res in js["results"]:
        t = res["thermo"]
        out[res["name"]] = dict(E=t["scf_energy"], Esp=t["sp_energy"], ZPE=t["zpe"], H=t["enthalpy"],
                                TS=TEMP * t["entropy"], TqhS=TEMP * t["qh_entropy"], G=t["gibbs_free_energy"],
                                qhG=t["qh_gibbs_free_energy"], imag=list(t.get("im_frequency_wn") or []),
                                nfreq=len(t.get("frequency_wn") or []))
    return out, " ".join(["goodvibes", *args])


# ------------------------------------------------------------------------------------------------ build
def main():
    gv_version = check_goodvibes()
    built = date.today().isoformat()
    species_cfg = list(C.SPECIES)
    for irc in C.IRCS:
        for e in irc["ends"]:
            species_cfg.append(dict(e, group="irc", irc=irc["id"], kind=e.get("kind", "minimum"),
                                    keys=e.get("keys", irc["keys"]),
                                    src=f"IRC_Section: Opt + Freq (TightOPT) from the {e['dir']} end of the "
                                        f"{irc['label'].replace('IRC: ', '')} IRC"))

    # --- read every structure
    sp, problems, method_rows = {}, [], []
    for s in species_cfg:
        d = os.path.join(CALC, s["id"])
        if not os.path.isfile(os.path.join(d, "output.out.gz")):
            continue
        text = read_gz(os.path.join(d, "output.out.gz"))
        normal = "ORCA TERMINATED NORMALLY" in text
        converged = "THE OPTIMIZATION HAS CONVERGED" in text
        freqs = last_freqs(text) if normal else None
        status = ("complete" if normal and freqs else "freq" if converged else "unconverged")
        ok, mc = method_check(text, "opt")
        method_rows.append(dict(id=s["id"], what="M06-2X Opt/OptTS + Freq", ok=ok, **mc, status=status))
        if not ok:
            problems.append(f"{s['id']}: method differs from the SI: {', '.join(mc['failed'])}")
        atoms = last_geometry(text)
        rec = dict(id=s["id"], label=s["label"], kind=s["kind"], group=s.get("group", "main"), status=status,
                   E=last_energy(text), atoms=atoms, natoms=len(atoms), keys=[list(k) for k in s.get("keys", [])],
                   source=s.get("src") or "", freq_type=mc["freq"], irc=s.get("irc"), dir=s.get("dir"),
                   match=s.get("match"), classify=s.get("classify"),
                   forming=s.get("forming", []), breaking=s.get("breaking", []))
        if os.path.isfile(os.path.join(d, "dlpno.out.gz")):
            ctext = read_gz(os.path.join(d, "dlpno.out.gz"))
            okc, mcc = method_check(ctext, "dlpno")
            method_rows.append(dict(id=s["id"], what="DLPNO-CCSD(T) single point", ok=okc, **mcc,
                                    status="complete" if "ORCA TERMINATED NORMALLY" in ctext else "running"))
            if not okc:
                problems.append(f"{s['id']} DLPNO: method differs from the SI: {', '.join(mcc['failed'])}")
            ci = dlpno_info(ctext)
            if ci["normal"]:
                rec["E_cc"], rec["T1"] = ci["E"], ci["T1"]
        sp[s["id"]] = rec

    # IRC runs
    irc_runs = {}
    for irc in C.IRCS:
        d = os.path.join(CALC, irc["id"])
        if not os.path.isfile(os.path.join(d, "output.out.gz")):
            continue
        text = read_gz(os.path.join(d, "output.out.gz"))
        ok, mc = method_check(text, "irc")
        method_rows.append(dict(id=irc["id"], what="M06-2X IRC", ok=ok, **mc,
                                status="complete" if "ORCA TERMINATED NORMALLY" in text else "running"))
        if not ok:
            problems.append(f"{irc['id']}: method differs from the SI: {', '.join(mc['failed'])}")
        irc_runs[irc["id"]] = dict(converged=text.count("THE IRC HAS CONVERGED"),
                                   normal="ORCA TERMINATED NORMALLY" in text)

    # --- GoodVibes: 1 mol/L, 1 atm, and composite with DLPNO-CCSD(T) where available
    with_freq = [k for k, v in sp.items() if v["status"] == "complete"]
    gv1, cmd1 = run_goodvibes(with_freq, "1M", conc=True)
    gva, cmda = run_goodvibes(with_freq, "1atm", conc=False)
    with_cc = [k for k in with_freq if "E_cc" in sp[k]]
    gvc, cmdc = run_goodvibes(with_cc, "DLPNO_1M", conc=True, spc=True) if with_cc else ({}, "")
    ss = None
    for k in with_freq:
        r1, ra = gv1[k], gva[k]
        assert abs(r1["E"] - sp[k]["E"]) < 2e-6, (k, r1["E"], sp[k]["E"])
        s = sp[k]
        s.update(ZPE=r1["ZPE"], H=r1["H"], TS=r1["TS"], TqhS=r1["TqhS"], G_rrho=r1["G"], G=r1["qhG"], G_atm=ra["qhG"],
                 imag=r1["imag"], nimag=len(r1["imag"]), Gcorr=r1["qhG"] - r1["E"])
        d_ss = (r1["qhG"] - ra["qhG"]) * KJ
        ss = d_ss if ss is None else ss
        assert abs(d_ss - ss) < 1e-3, ("standard-state correction differs", k, d_ss, ss)
        if k in gvc:
            s["Gc"] = gvc[k]["qhG"]
            assert abs(gvc[k]["Esp"] - s["E_cc"]) < 2e-6, k
            assert abs(s["Gc"] - (s["E_cc"] + s["Gcorr"])) < 2e-6, (k, "composite G")
        want = 1 if s["kind"] == "ts" else 0
        if s["nimag"] != want:
            problems.append(f"{k}: {s['nimag']} imaginary frequencies ({s['imag']}), expected {want}")
    expected_ss = 8.314462618 * TEMP / 1000 * math.log(8.314462618 * TEMP / 101325 * 1000)
    assert abs(ss - expected_ss) < 0.01, (ss, expected_ss)

    # --- IRC end points: which side a CO-IRC end is, and comparison with the structure it should be
    alias = {}
    for k, s in sp.items():
        if s.get("classify") == "co":
            dco = (dist(s["atoms"], 21, 28) + dist(s["atoms"], 25, 28)) / 2
            role = "prodCO" if dco > 2.3 else "adduct"
            s["match"] = "complex_prodCO" if role == "prodCO" else "adduct"
            s["label"] = ("Decarbonylated product·CO complex (CO-extrusion IRC end)" if role == "prodCO"
                          else "[4+2] adduct (CO-extrusion IRC end)")
            s["kind"] = "complex" if role == "prodCO" else "minimum"
            alias["ep_CO_" + role] = k
    for k, s in sp.items():
        m = s.get("match")
        if s["group"] != "irc" or not m or m not in sp:
            continue
        measure, val, heavy = compare(s["atoms"], sp[m]["atoms"])
        dE = (s["E"] - sp[m]["E"]) * KJ
        dG = (s["G"] - sp[m]["G"]) * KJ if s.get("G") is not None and sp[m].get("G") is not None else None
        thresh = 0.3 if measure == "RMSD" else 0.03
        same = val is not None and val < thresh and abs(dE) < 0.5
        s["compare"] = dict(match=m, measure=measure, value=val, heavy=heavy, dE=dE, dG=dG, same=same)

    def resolve(c):
        return alias.get(c, c)

    def pick(cands):
        avail = [resolve(c) for c in cands if resolve(c) in sp]
        for c in avail:
            if sp[c].get("G") is not None:
                return c
        for c in avail:
            if sp[c].get("E") is not None:
                return c
        return None

    # --- states (relative to carbene + tetracyclone), reactions (relative to their first row)
    def total(ids, key):
        vals = [sp[i].get(key) for i in ids]
        return None if any(v is None for v in vals) else sum(vals)

    states = {}
    for sid, (label, members, short) in C.STATES.items():
        used = [pick(c) for c in members]
        if None in used:
            continue
        fallback = [c[0] for c, u in zip(members, used) if u != resolve(c[0])]
        states[sid] = dict(id=sid, label=label, short=short, members=used, fallback_for=fallback,
                           kind="ts" if any(sp[m]["kind"] == "ts" for m in used) else
                                ("complex" if any(sp[m]["kind"] == "complex" for m in used) else "minimum"),
                           abs=dict(E=total(used, "E"), G=total(used, "G"), G_atm=total(used, "G_atm"),
                                    Gc=total(used, "Gc"), E_cc=total(used, "E_cc")))
    ref = states["R"]["abs"]

    def rel(a, b, key):
        return None if a[key] is None or b[key] is None else round((a[key] - b[key]) * KJ, 2)

    def formula(ids):
        f = {}
        for i in ids:
            for a in sp[i]["atoms"]:
                f[a[0]] = f.get(a[0], 0) + 1
        return f

    f_ref = formula(states["R"]["members"])
    for st in states.values():   # on the PDF's scale only if the state has the same atoms as carbene + tetracyclone
        same = formula(st["members"]) == f_ref
        st.update(E=rel(st["abs"], ref, "E") if same else None, G=rel(st["abs"], ref, "G") if same else None,
                  G_atm=rel(st["abs"], ref, "G_atm") if same else None, Gc=rel(st["abs"], ref, "Gc") if same else None)
    reactions = []
    for rx in C.REACTIONS:
        rows = [r for r in rx["rows"] if r in states]
        if not rows:
            continue
        assert all(formula(states[r]["members"]) == formula(states[rows[0]]["members"]) for r in rows), rx["id"]
        r0 = states[rows[0]]["abs"]
        reactions.append(dict(id=rx["id"], label=rx["label"], rows=[
            dict(state=r, label=states[r]["label"], kind=states[r]["kind"], members=states[r]["members"],
                 dE=rel(states[r]["abs"], r0, "E"), dG=rel(states[r]["abs"], r0, "G"),
                 dG_atm=rel(states[r]["abs"], r0, "G_atm"), dGc=rel(states[r]["abs"], r0, "Gc"),
                 dEcc=rel(states[r]["abs"], r0, "E_cc")) for r in rows]))

    # --- TS table and IRC table
    ts_table = []
    for k, s in sp.items():
        if s["kind"] != "ts":
            continue
        ts_table.append(dict(id=k, label=s["label"], imag=s.get("imag", []),
                             forming=[dict(bond=f"C{i}–C{j}", d=round(dist(s["atoms"], i, j), 3)) for i, j in s["forming"]],
                             breaking=[dict(bond=f"C{i}–C{j}", d=round(dist(s["atoms"], i, j), 3)) for i, j in s["breaking"]]))
    ircs = []
    for irc in C.IRCS:
        ends = []
        for e in irc["ends"]:
            s = sp.get(e["id"])
            ends.append(dict(id=e["id"], dir=e["dir"], present=s is not None,
                             label=s["label"] if s else e["label"], status=s["status"] if s else "not started",
                             nimag=s.get("nimag") if s else None, compare=s.get("compare") if s else None,
                             dE_ts=round((s["E"] - sp[irc["ts"]]["E"]) * KJ, 1) if s else None))
        run = irc_runs.get(irc["id"])
        ircs.append(dict(id=irc["id"], label=irc["label"], ts=irc["ts"], done=run is not None,
                         converged=run["converged"] if run else 0, ends=ends))

    # --- paths
    def path_ref(r):
        return sp["carbene"]["E"] if r == "carbene" else ref["E"]
    paths = []
    plist = [dict(id=i["id"], label=i["label"], group="IRCs (IRC_Section)", keys=i["keys"], ref=i["ref"],
                  file=os.path.join(CALC, i["id"], "trajectory.xyz"),
                  note=f"IRC in both directions from {sp[i['ts']]['label']}; each end was then optimised with "
                       f"frequencies (see the Energies tab).") for i in C.IRCS]
    plist += [dict(p, file=os.path.join(CALC, "paths", p["id"], "trajectory.xyz")) for p in C.PATHS]
    for p in plist:
        if not os.path.isfile(p["file"]):
            continue
        frames = read_frames(p["file"])
        r = path_ref(p["ref"])
        energies = [round((comment_energy(c) - r) * KJ, 2) if comment_energy(c) is not None else None for c, _ in frames]
        paths.append(dict(id=p["id"], label=p["label"], group=p["group"], keys=[list(k) for k in p["keys"]],
                          note=p["note"], natoms=len(frames[0][1]), energies=energies, file=f"xyz/paths/{p['id']}.xyz",
                          ref="the [2.2.2] carbene" if p["ref"] == "carbene" else "carbene + tetracyclone",
                          frames=[xyz_text(a, f"{p['label']} frame {i + 1}") for i, (c, a) in enumerate(frames)]))

    # --- pending work, as shown on the site
    pending = []
    for i in ircs:
        if not i["done"]:
            pending.append(f"{i['label']}: the IRC is still running; Opt + Freq of its two ends follow")
        for e in i["ends"]:
            if e["present"] and e["status"] == "freq":
                pending.append(f"{e['label']}: optimised; frequencies still running")
            elif e["present"] and e["status"] == "unconverged":
                pending.append(f"{e['label']}: optimisation still running")
    no_g = [st["label"] for sid, st in states.items() if st["abs"]["G"] is None and sid in
            {s for r in C.ROUTES for s in r["states"]}]
    if no_g:
        pending.append("No free energy yet (electronic energy only) for: " + "; ".join(no_g))
    n_cc_needed = len({m for st in states.values() for m in st["members"]})
    n_cc = len({m for st in states.values() for m in st["members"] if "E_cc" in sp[m]})
    pending.append(f"DLPNO-CCSD(T) single points: {n_cc} of the {n_cc_needed} structures in the diagrams are done; "
                   "relative energies stay at the M06-2X level until all are")

    # --- write everything
    for d in ("data", "xyz/species", "xyz/paths", "thermochemistry"):
        os.makedirs(os.path.join(REPO, d), exist_ok=True)
    for f in os.listdir(os.path.join(REPO, "xyz/species")):   # regenerate from scratch
        os.remove(os.path.join(REPO, "xyz/species", f))
    for f in os.listdir(os.path.join(REPO, "xyz/paths")):
        os.remove(os.path.join(REPO, "xyz/paths", f))
    out_species = []
    for k, s in sp.items():
        g = s.get("G")
        head = (f"{s['label']} | M06-2X/def2-TZVP D3(0) CPCM(CH2Cl2) | Eel(M06-2X) = {s['E']:.6f} Eh"
                + (f"; Eel(DLPNO-CCSD(T)) = {s['E_cc']:.6f} Eh" if "E_cc" in s else "")
                + (f"; G(M06-2X, 313.15 K, 1 M) = {g:.6f} Eh" if g is not None else "; G: frequencies pending"))
        fname = f"xyz/species/{k}.xyz"
        with open(os.path.join(REPO, fname), "w") as fh:
            fh.write(xyz_text(s["atoms"], head))
        o = dict(id=k, label=s["label"], kind=s["kind"], group=s["group"], status=s["status"], natoms=s["natoms"],
                 E_Eh=round(s["E"], 9), G_Eh=round(g, 6) if g is not None else None,
                 G_atm_Eh=round(s["G_atm"], 6) if s.get("G_atm") is not None else None,
                 Gcorr_Eh=round(s["Gcorr"], 6) if s.get("Gcorr") is not None else None,
                 ZPE_Eh=round(s["ZPE"], 6) if s.get("ZPE") is not None else None,
                 H_Eh=round(s["H"], 6) if s.get("H") is not None else None,
                 TS_Eh=round(s["TS"], 6) if s.get("TS") is not None else None,
                 TqhS_Eh=round(s["TqhS"], 6) if s.get("TqhS") is not None else None,
                 Ecc_Eh=round(s["E_cc"], 6) if "E_cc" in s else None, T1=round(s["T1"], 3) if s.get("T1") else None,
                 Gc_Eh=round(s["Gc"], 6) if "Gc" in s else None,
                 imag=[round(x, 2) for x in s.get("imag", [])], nimag=s.get("nimag"), freq_type=s["freq_type"],
                 keys=s["keys"], source=s["source"], irc=s["irc"], dir=s["dir"], compare=s.get("compare"),
                 file=fname, xyz=xyz_text(s["atoms"], head),
                 states=[sid for sid, st in states.items() if k in st["members"]])
        out_species.append(o)
    for p in paths:
        with open(os.path.join(REPO, p["file"]), "w") as fh:
            for i, fr in enumerate(p["frames"]):
                e = p["energies"][i]
                lines = fr.split("\n")
                lines[1] = f"{p['label']} frame {i + 1}" + (f" | {e:+.2f} kJ/mol rel. to {p['ref']}" if e is not None else "")
                fh.write("\n".join(lines))
    meta = dict(
        title="FBW tetracyclone trapping explorer", built=built, url=SITE_URL,
        level="M06-2X/def2-TZVP D3(0), CPCM(CH₂Cl₂: ε = 9.08, n = 1.424), RIJCOSX (def2/J), DefGrid2, TightSCF, "
              "TightOPT; ORCA 6.1.0",
        thermo=f"GoodVibes {gv_version} at 313.15 K from the M06-2X harmonic frequencies: Grimme quasi-RRHO entropy "
               "(100 cm⁻¹ cut-off), no quasi-harmonic enthalpy correction, no frequency scaling; 1 mol/L standard state",
        gv_command=cmd1 + " *.out", gv_command_atm=cmda + " *.out", gv_command_cc=(cmdc + " *.out") if cmdc else "",
        ss=round(ss, 2), eh=KJ, temp=TEMP, gv_version=gv_version,
        cc_level="DLPNO-CCSD(T)/def2-TZVPP (def2-TZVPP/C, TightPNO), same CPCM, on the M06-2X geometries",
        reference="All energies in kJ/mol relative to the [2.2.2] carbene and tetracyclone as separate molecules",
        cc_check={r["label"]: r["dGc"] for rx in reactions if rx["id"] == "fbw" for r in rx["rows"][1:]
                  if r["dGc"] is not None},
        n_cc=n_cc, n_cc_needed=n_cc_needed, pending=pending, problems=problems)
    with open(os.path.join(REPO, "data/data.js"), "w") as fh:
        fh.write("// Generated by tools/build_site.py; do not edit by hand.\n")
        fh.write("window.SB_DATA = " + json.dumps(dict(meta=meta, species=out_species, states=states, routes=C.ROUTES,
                                                        reactions=reactions, ts_table=ts_table, ircs=ircs,
                                                        methods=method_rows), ensure_ascii=False, indent=1) + ";\n")
    with open(os.path.join(REPO, "data/paths.js"), "w") as fh:
        fh.write("// Generated by tools/build_site.py; do not edit by hand.\n")
        fh.write("window.SB_PATHS = " + json.dumps(paths, ensure_ascii=False) + ";\n")
    write_tables(out_species, reactions, ts_table, ircs, method_rows, meta)
    write_readme(out_species, states, reactions, ts_table, ircs, paths, meta)
    open(os.path.join(REPO, ".nojekyll"), "w").close()
    print(f"structures {len(sp)} (with frequencies {len(with_freq)}, with DLPNO-CCSD(T) {len(with_cc)}); "
          f"states {len(states)}; paths {len(paths)}; standard-state correction {ss:.3f} kJ/mol")
    for p in pending:
        print("  pending:", p)
    for p in problems:
        print("  PROBLEM:", p)


# ------------------------------------------------------------------------------------------------ tables + README
def f6(v):
    return "" if v is None else f"{v:.6f}"


def f1(v, signed=True):
    if v is None:
        return "–"
    return f"{v:+.1f}" if signed else f"{v:.1f}"


def bonds(lst):
    return ", ".join(f"{b['d']:.3f}" for b in lst) or "—"


def write_tables(species, reactions, ts_table, ircs, methods, meta):
    T = THERMO
    with open(os.path.join(T, "energies_absolute.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["id", "structure", "type", "status", "Nimag", "imag_freq_cm-1", "Eel_M06-2X_Eh", "ZPE_Eh", "H_Eh",
                    "T.S_Eh", "T.qh-S_Eh", "G_minus_Eel_Eh", "G_M06-2X_1M_Eh", "G_M06-2X_1atm_Eh",
                    "Eel_DLPNO-CCSD(T)_Eh", "T1", "G_composite_1M_Eh"])
        for s in species:
            w.writerow([s["id"], s["label"], s["kind"], s["status"], "" if s["nimag"] is None else s["nimag"],
                        " ".join(f"{x:.2f}" for x in s["imag"]), f"{s['E_Eh']:.6f}", f6(s["ZPE_Eh"]), f6(s["H_Eh"]),
                        f6(s["TS_Eh"]), f6(s["TqhS_Eh"]), f6(s["Gcorr_Eh"]), f6(s["G_Eh"]), f6(s["G_atm_Eh"]),
                        f6(s["Ecc_Eh"]), "" if s["T1"] is None else f"{s['T1']:.3f}", f6(s["Gc_Eh"])])
    with open(os.path.join(T, "energies_relative.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["reaction", "state", "structures", "dEel_kJ/mol", "dG_1M_kJ/mol", "dG_1atm_kJ/mol",
                    "dEel_DLPNO_kJ/mol", "dG_composite_1M_kJ/mol"])
        for rx in reactions:
            for r in rx["rows"]:
                w.writerow([rx["label"], r["label"], " + ".join(r["members"]), r["dE"], r["dG"], r["dG_atm"],
                            r["dEcc"], r["dGc"]])
    with open(os.path.join(T, "transition_structures.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["id", "transition structure", "imag_freq_cm-1", "forming", "forming_A", "breaking", "breaking_A"])
        for t in ts_table:
            w.writerow([t["id"], t["label"], " ".join(f"{x:.2f}" for x in t["imag"]),
                        " ".join(b["bond"] for b in t["forming"]), bonds(t["forming"]),
                        " ".join(b["bond"] for b in t["breaking"]), bonds(t["breaking"])])
    with open(os.path.join(T, "irc_end_points.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["irc", "TS", "IRC done", "end", "direction", "status", "Nimag", "E_rel_TS_kJ/mol", "compared_with",
                    "measure", "value_A", "heavy_atom_RMSD_A", "dEel_kJ/mol", "dG_kJ/mol", "same_structure"])
        for i in ircs:
            for e in i["ends"]:
                c = e["compare"] or {}
                w.writerow([i["id"], i["ts"], i["done"], e["id"], e["dir"], e["status"],
                            "" if e["nimag"] is None else e["nimag"], "" if e["dE_ts"] is None else e["dE_ts"],
                            c.get("match", ""), c.get("measure", ""),
                            "" if c.get("value") is None else f"{c['value']:.3f}",
                            "" if c.get("heavy") is None else f"{c['heavy']:.3f}",
                            "" if c.get("dE") is None else f"{c['dE']:.2f}",
                            "" if c.get("dG") is None else f"{c['dG']:.2f}", c.get("same", "")])
    with open(os.path.join(T, "method_check.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["id", "calculation", "ORCA version", "method line(s)", "frequencies", "matches the SI",
                    "differences", "status"])
        for m in methods:
            w.writerow([m["id"], m["what"], m["version"], m["method_line"], m["freq"], m["ok"], "; ".join(m["failed"]),
                        m["status"]])
    # Cartesian coordinates in the SI's format
    with open(os.path.join(T, "SI_coordinates.txt"), "w") as fh:
        fh.write("Cartesian coordinates (Å) at the M06-2X/def2-TZVP/D3(0)/CPCM(CH2Cl2) level. All species are closed-shell "
                 "singlets (charge 0, multiplicity 1). G is the free energy at 313.15 K and 1 mol/L (M06-2X until the "
                 "DLPNO-CCSD(T) energies are complete).\n\n")
        for s in species:
            parts = [f"Eel(M06-2X) = {s['E_Eh']:.6f} Eh"]
            if s["Ecc_Eh"] is not None:
                parts.append(f"Eel(DLPNO-CCSD(T)) = {s['Ecc_Eh']:.6f} Eh")
            parts.append(f"G = {s['G_Eh']:.6f} Eh" if s["G_Eh"] is not None else "G: frequencies pending")
            fh.write(f"{s['label']}\n" + ";  ".join(parts) + "\n")
            fh.write("\n".join(s["xyz"].split("\n")[2:]).rstrip() + "\n\n")


def write_readme(species, states, reactions, ts_table, ircs, paths, meta):
    def rx_table(rx):
        rows = "\n".join(f"| {r['label']} | {f1(r['dE'])} | {f1(r['dG'])} | {f1(r['dG_atm'])} |" for r in rx["rows"])
        return f"**{rx['label']}**\n\n| Structure | ΔEel | ΔG (1 M) | ΔG (1 atm) |\n| --- | --- | --- | --- |\n{rows}\n"
    rx_md = "\n".join(rx_table(rx) for rx in reactions if rx["id"] != "all")
    all_rx = next(rx for rx in reactions if rx["id"] == "all")
    ts_md = "\n".join(f"| {t['label']} | {', '.join(f'{x:.1f}' for x in t['imag']) or '–'} | {bonds(t['forming'])} | "
                      f"{bonds(t['breaking'])} |" for t in ts_table)
    irc_rows = []
    for i in ircs:
        for e in i["ends"]:
            c = e["compare"]
            cmp = "–"
            if c:
                cmp = (f"{'same as' if c['same'] else 'differs from'} `{c['match']}` "
                       f"({c['measure']} {c['value']:.3f} Å, ΔE {c['dE']:+.2f} kJ/mol"
                       + (f", ΔG {c['dG']:+.2f}" if c["dG"] is not None else "") + ")")
            irc_rows.append(f"| {i['label']} | {e['dir']} | {e['label'] if e['present'] else '–'} | "
                            f"{e['status']} | {'' if e['nimag'] is None else e['nimag']} | {cmp} |")
    sp_md = "\n".join(
        f"| {s['label']} | {s['kind']} | {'' if s['nimag'] is None else s['nimag']} | {s['E_Eh']:.6f} | "
        f"{f6(s['Gcorr_Eh']) or '–'} | {f6(s['G_Eh']) or '–'} | {f6(s['Ecc_Eh']) or '–'} | "
        f"{'' if s['T1'] is None else s['T1']} | `{s['file']}` |" for s in species)
    pending = "\n".join(f"- {p}" for p in meta["pending"])
    problems = ("\n**Checks that failed:**\n\n" + "\n".join(f"- {p}" for p in meta["problems"]) + "\n") if meta["problems"] else ""
    text = f"""# FBW tetracyclone trapping explorer

An interactive viewer for the reactants, intermediates, transition states and products of the 6→7
Fritsch–Buttenberg–Wiechell (FBW) system trapped by tetracyclone (section B of the study's ORCA work summary).
The [2.2.2] carbene either adds directly to tetracyclone ([1+2], two faces), or first rearranges to the [3.2.2] alkyne
(FBW), which tetracyclone traps by a Diels–Alder reaction followed by CO extrusion.

**Live site:** {meta['url']}. You can also open `index.html` directly in a browser: no server or installation is
needed, and the page makes no network connections.

## What you can do

- **Routes:** energy diagrams of each route (free energy at 1 mol/L, or electronic energy). Click a level to see its structure in 3D, with
  the bonds that form or break shown as dashed lines with their lengths.
- **Structures:** every structure, including the optimised IRC end points, with its energies and imaginary frequency.
- **Paths:** step through the IRCs, NEB bands and relaxed scans frame by frame, with the energy along the path.
- **Energies:** the tables of the Supplementary Information for this section: relative energies by reaction, the
  transition structures, the IRC connections, the absolute energies and a check of every calculation's method.
- **Download:** every structure is in `xyz/`, and the tables are in `thermochemistry/`.

## Method (as in Part II of the Supplementary Information)

- **Structures:** {meta['level']}. Minima have no imaginary frequency and transition structures exactly one.
  `thermochemistry/method_check.csv` compares every input with this.
- **Thermochemistry:** {meta['thermo']}. The 1 mol/L values add RT ln(Vm/V°) = {meta['ss']:.2f} kJ/mol to the
  1 atm value of each species. Commands (run by `tools/build_site.py`):
  - `{meta['gv_command']}`
  - `{meta['gv_command_atm']}` (1 atm values)
- **Single points:** {meta['cc_level']}. Composite G = Eel(DLPNO-CCSD(T)) + [G − Eel](M06-2X), from GoodVibes `--spc`.
  {meta['n_cc']} of the {meta['n_cc_needed']} structures in the diagrams have their DLPNO-CCSD(T) energy so far, so
  **all relative energies below are M06-2X values** for now.
- **Units:** kJ/mol, 1 Eh = {meta['eh']:.2f} kJ/mol, relative to the separated reactants of each step.
{problems}
## Relative energies by reaction (kJ/mol, M06-2X, 313.15 K)

{rx_md}
## Whole surface, relative to [2.2.2] carbene + tetracyclone (the PDF's 0)

| State | ΔEel | ΔG (1 M) | ΔG (1 atm) | Structures used |
| --- | --- | --- | --- | --- |
{chr(10).join(f"| {r['label']} | {f1(r['dE'])} | {f1(r['dG'])} | {f1(r['dG_atm'])} | {' + '.join(f'`{m}`' for m in r['members'])} |" for r in all_rx['rows'])}

– means no frequencies yet. An IRC end point is used as soon as its frequencies are done; until then the earlier
structure (or the end point's electronic energy alone) stands in.

## Transition structures

| Transition structure | ν_imag (cm⁻¹) | Forming (Å) | Breaking (Å) |
| --- | --- | --- | --- |
{ts_md}

FBW shift: distance from the migrating carbon (C5) to its destination (C16, forming) and origin (C2, breaking).
Additions: the two C···C distances between the reacting fragments. CO extrusion: the two breaking C–C(O) bonds.

## IRCs and their optimised end points

| IRC | Direction | End point | State | Nimag | Compared with |
| --- | --- | --- | --- | --- | --- |
{chr(10).join(irc_rows)}

## Still to come

{pending}

## Absolute energies (Eh)

| Structure | Type | Nimag | Eel M06-2X | G − Eel | G (M06-2X, 1 M) | Eel DLPNO-CCSD(T) | T1 | File |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
{sp_md}

## Files

- `index.html`, `assets/`: the viewer (3Dmol.js, BSD-3-Clause; see `assets/3Dmol_LICENSE.txt`).
- `data/`: the structures, energies and paths as script files for the viewer.
- `xyz/species/`, `xyz/paths/`: the structures and paths as xyz files.
- `thermochemistry/`: the GoodVibes output (`GoodVibes_1M.dat`, `GoodVibes_1atm.dat`, `GoodVibes_DLPNO_1M.dat`, and
  JSON), the tables above as CSV, and `SI_coordinates.txt` (coordinates with energies, in the SI's format).
- `calculations/<id>/`: the ORCA outputs everything is computed from (`output.out.gz`, `dlpno.out.gz`), with each
  input and where it came from. Host names and working directories are redacted.
- `tools/`: `site_config.py` (what is shown), `collect.py` (copies new results from the project folder) and
  `build_site.py` (runs GoodVibes and rebuilds the site).

## Rebuilding

```
pip install "goodvibes==4.3.0"
python3 tools/collect.py --root HOST:/path/to/project   # only to fetch new results
python3 tools/build_site.py
```

Built {meta['built']}.
"""
    open(os.path.join(REPO, "README.md"), "w").write(text)


if __name__ == "__main__":
    main()
