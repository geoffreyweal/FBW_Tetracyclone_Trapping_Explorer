"""What the explorer shows: structures, IRCs, paths, energy-diagram states, routes and reaction tables.

Source paths are relative to the project folder on the cluster (the folder holding attempt1/, attempt2/, attempt3/,
IRC_Section/, DLPNO_Section/ and original_data/). tools/collect.py copies them into calculations/<id>/, and
tools/build_site.py reads only calculations/, so the site can be rebuilt from this repository alone.

Atom indices are ORCA's, counting from 0, and are the same in every 71-atom structure of a route.
"""

W = "original_data/FBW_paper/working"
F = W + "/06_analysis/final_DFT"              # the paper's Opt+Freq outputs, as used for its GoodVibes table
V2 = "attempt1/orca6_TSopt_v2"                 # attempt-1 TSs (OptTS + Freq, analytic Hessian)
S1 = "attempt3/Step1_minima"                   # attempt-3 pre-reaction complexes
LN = "attempt3/Step3_TS_from_long-NEB"         # attempt-3 TS searches from the long-NEB maxima
IRC = "IRC_Section/Step1_IRC"                  # IRCs from the five TSs
EP = "IRC_Section/Step2_IRC_endpoints"         # Opt + Freq of both ends of each IRC
CC1 = "DLPNO_Section/Step1_DLPNO"              # DLPNO-CCSD(T) single points: TSs and separate molecules
CC2 = "DLPNO_Section/Step2_DLPNO_IRC_endpoints"  # DLPNO-CCSD(T) single points: optimised IRC end points

# Bonds that form or break (shown dashed in 3D, and used in the TS table)
K_Z = [(15, 22), (15, 21)]      # [1+2], Z face: C15 (carbene C) adds to C22 first
K_E = [(15, 21), (15, 22)]      # [1+2], E face: C15 adds to C21 first
K_FBW = [(2, 5), (5, 16)]       # FBW: C5 migrates from C2 (breaking) to C16 (forming)
K_DA = [(14, 21), (18, 25)]     # [4+2]: alkyne C14/C18 add to diene termini C21/C25
K_CO = [(21, 28), (25, 28)]     # CO extrusion: C28(=O) leaves C21 and C25

# ------------------------------------------------------------------------------------------------ structures
# dft: the M06-2X output (Opt/OptTS + Freq) giving the energy, geometry and frequencies for GoodVibes.
# cc:  DLPNO-CCSD(T) single point(s) on that geometry; the first finished one is used.
# forming/breaking: bonds for the transition-structure table (SI definitions: for the FBW shift, the migrating
#   carbon's distance to its destination and to its origin; for additions, the two C···C distances between the
#   reacting fragments).
SPECIES = [
    dict(id="carbene", label="[2.2.2] carbene", kind="minimum", dft=f"{F}/67-carbene.out",
         cc=[f"{CC1}/67-carbene/output.out", f"{F}/67-carbene_CC.out"], src="Original study (DFT minimum)"),
    dict(id="tetracyclone", label="Tetracyclone (trapping reagent)", kind="minimum", dft=f"{F}/tetracyclone.out",
         cc=[f"{CC1}/tetracyclone/output.out", "attempt1/orca6_DLPNO/tetracyclone/output.out"],
         src="Original study (DFT minimum)"),
    dict(id="complex_Z", label="Carbene·tetracyclone complex, Z face", kind="complex",
         dft=f"{S1}/path1a_carbene-tetracyclone-complex-Z/output.out", keys=K_Z,
         src="Attempt 3, Step 1 (Opt + Freq), for the NEB runs"),
    dict(id="TS_SSZ", label="S,S,Z-#-TSa ([1+2], Z face)", kind="ts",
         dft=f"{V2}/67-carbene-tetracyclone-SSZ-trap-TSa/output.out",
         cc=[f"{CC1}/67-carbene-tetracyclone-SSZ-trap-TSa/output.out"], keys=K_Z, forming=K_Z,
         src="Attempt 1 (OptTS + Freq, analytic Hessian, v2); reproduced by the attempt-3 NEB-TS and the attempt-2 scan"),
    dict(id="SSZ", label="S,S,Z-# ([1+2] adduct)", kind="minimum", dft=f"{F}/67-carbene-tetracyclone-SSZ-trap.out",
         keys=K_Z, src="Original study (DFT minimum)"),
    dict(id="complex_E", label="Carbene·tetracyclone complex, E face", kind="complex",
         dft=f"{S1}/path1b_carbene-tetracyclone-complex-E/output.out", keys=K_E,
         src="Attempt 3, Step 1 (Opt + Freq), for the NEB runs"),
    dict(id="TS_SSE", label="S,S,E-#-TSb ([1+2], E face)", kind="ts",
         dft=f"{V2}/67-carbene-tetracyclone-SSE-trap-TSb/output.out",
         cc=[f"{CC1}/67-carbene-tetracyclone-SSE-trap-TSb/output.out"], keys=K_E, forming=K_E,
         src="Attempt 1 (OptTS + Freq, analytic Hessian, v2); reproduced by the attempt-3 NEB-TS and the attempt-2 scan"),
    dict(id="SSE", label="S,S,E-# ([1+2] adduct)", kind="minimum", dft=f"{F}/67-carbene-tetracyclone-SSE-trap.out",
         keys=K_E, src="Original study (DFT minimum)"),
    dict(id="TS_FBW", label="#-FBW-TS (FBW rearrangement)", kind="ts", dft=f"{F}/67-FBW-TS.out",
         cc=[f"{CC1}/67-FBW-TS/output.out", f"{F}/67-FBW-TS_CC.out"], keys=K_FBW, forming=[(5, 16)], breaking=[(2, 5)],
         src="Original study (OptTS + numerical frequencies); reproduced by the attempt-3 NEB-TS (0.0002 Å)"),
    dict(id="alkyne", label="[3.2.2] alkyne", kind="minimum", dft=f"{F}/67-alkyne.out",
         cc=[f"{CC1}/67-alkyne/output.out", f"{F}/67-alkyne_CC.out"], src="Original study (DFT minimum)"),
    dict(id="complex_alkyne", label="Alkyne·tetracyclone complex (attempt 3)", kind="complex",
         dft=f"{S1}/path2b_alkyne-tetracyclone-complex/output.out", keys=K_DA,
         src="Attempt 3, Step 1 (optimisation converged; cancelled during the frequencies). The [4+2] IRC end point "
             "is the same structure and replaces it"),
    dict(id="TS_DA", label="#-TS [4+2] (Diels–Alder)", kind="ts", dft=f"{LN}/path2b_DielsAlder-TS/output.out",
         cc=[f"{CC1}/67-alkyne-tetracyclone-trap-TS/output.out"], keys=K_DA, forming=K_DA,
         src="Attempt 3 (OptTS + Freq, analytic Hessian, from the [4+2] maximum of the long NEB); the lowest [4+2] TS "
             "found, with the alkyne cage turned differently from the attempt-1 TS"),
    dict(id="TS_DA_a1", label="#-TS [4+2], higher conformer (attempt 1)", kind="ts",
         dft=f"{V2}/67-alkyne-tetracyclone-trap-TS/output.out", keys=K_DA, forming=K_DA,
         src="Attempt 1 (OptTS + Freq, analytic Hessian, v2); matched by the attempt-2 scan maximum. Not used in the "
             "diagrams"),
    dict(id="adduct", label="[4+2] adduct", kind="minimum", dft=f"{F}/67-alkyne-tetracyclone-trap.out", keys=K_CO,
         src="Original study (DFT minimum)"),
    dict(id="TS_CO", label="#-TS CO extrusion", kind="ts", dft=f"{V2}/67-alkyne-tetracyclone-decarbonylation-TS/output.out",
         cc=[f"{CC1}/67-alkyne-tetracyclone-decarbonylation-TS/output.out"], keys=K_CO, breaking=K_CO,
         src="Attempt 1 (OptTS + Freq, analytic Hessian, v2); matched by the attempt-2 scan maximum"),
    dict(id="complex_prodCO", label="Decarbonylated product·CO complex (attempt 3)", kind="complex",
         dft=f"{S1}/path2c_decarbonylated-product-CO-complex/output.out", keys=K_CO,
         src="Attempt 3, Step 1 (cancelled at cycle 53, nearly converged; no frequencies). To be replaced by the "
             "CO-extrusion IRC end point"),
    dict(id="product", label="Decarbonylated product", kind="minimum", dft=f"{F}/67-alkyne-tetracyclone-decarbonylated.out",
         cc=[f"{CC1}/67-alkyne-tetracyclone-decarbonylated/output.out"], src="Original study (DFT minimum)"),
    dict(id="CO", label="Carbon monoxide", kind="minimum", dft=f"{F}/carbon-monoxide.out",
         cc=[f"{CC1}/carbon-monoxide/output.out"], src="Original study (DFT minimum)"),
]

# ------------------------------------------------------------------------------------------------ IRCs
# Each IRC (both directions from its TS) and the Opt + Freq of each end. "match" is the structure the end is
# compared with; "classify" lets the builder decide which side an end is on from its geometry.
IRCS = [
    dict(id="irc_SSZ", label="IRC: S,S,Z [1+2]", ts="TS_SSZ", src=f"{IRC}/path1a_SSZ-TSa", keys=K_Z, ref="R",
         ends=[dict(id="ep_SSZ_F", dir="forward", src=f"{EP}/path1a_SSZ-TSa_forward",
                    cc=[f"{CC2}/path1a_SSZ-TSa_forward/output.out"], kind="complex", match="complex_Z",
                    label="Carbene·tetracyclone complex, Z face (S,S,Z IRC end)"),
               dict(id="ep_SSZ_B", dir="backward", src=f"{EP}/path1a_SSZ-TSa_backward",
                    cc=[f"{CC2}/path1a_SSZ-TSa_backward/output.out"], kind="minimum", match="SSZ",
                    label="S,S,Z-# (S,S,Z IRC end)")]),
    dict(id="irc_SSE", label="IRC: S,S,E [1+2]", ts="TS_SSE", src=f"{IRC}/path1b_SSE-TSb", keys=K_E, ref="R",
         ends=[dict(id="ep_SSE_F", dir="forward", src=f"{EP}/path1b_SSE-TSb_forward",
                    cc=[f"{CC2}/path1b_SSE-TSb_forward/output.out"], kind="minimum", match="SSE",
                    label="S,S,E-# (S,S,E IRC end)"),
               dict(id="ep_SSE_B", dir="backward", src=f"{EP}/path1b_SSE-TSb_backward",
                    cc=[f"{CC2}/path1b_SSE-TSb_backward/output.out"], kind="complex", match="complex_E",
                    label="Carbene·tetracyclone complex, E face (S,S,E IRC end)")]),
    dict(id="irc_FBW", label="IRC: FBW rearrangement", ts="TS_FBW", src=f"{IRC}/path2a_FBW-TS", keys=K_FBW, ref="carbene",
         ends=[dict(id="ep_FBW_F", dir="forward", src=f"{EP}/path2a_FBW-TS_forward",
                    cc=[f"{CC2}/path2a_FBW-TS_forward/output.out"], kind="minimum", match="carbene",
                    label="[2.2.2] carbene (FBW IRC end)"),
               dict(id="ep_FBW_B", dir="backward", src=f"{EP}/path2a_FBW-TS_backward",
                    cc=[f"{CC2}/path2a_FBW-TS_backward/output.out"], kind="minimum", match="alkyne",
                    label="[3.2.2] alkyne (FBW IRC end)")]),
    dict(id="irc_DA", label="IRC: [4+2] (lowest TS)", ts="TS_DA", src=f"{IRC}/path2b_DielsAlder-TS", keys=K_DA, ref="R",
         ends=[dict(id="ep_DA_F", dir="forward", src=f"{EP}/path2b_DielsAlder-TS_forward",
                    cc=[f"{CC2}/path2b_DielsAlder-TS_forward/output.out"], kind="minimum", match="adduct", keys=K_CO,
                    label="[4+2] adduct from the lowest [4+2] TS ([4+2] IRC end)"),
               dict(id="ep_DA_B", dir="backward", src=f"{EP}/path2b_DielsAlder-TS_backward",
                    cc=[f"{CC2}/path2b_DielsAlder-TS_backward/output.out"], kind="complex", match="complex_alkyne",
                    label="Alkyne·tetracyclone complex ([4+2] IRC end)")]),
    dict(id="irc_CO", label="IRC: CO extrusion", ts="TS_CO", src=f"{IRC}/path2c_CO-extrusion-TS", keys=K_CO, ref="R",
         ends=[dict(id="ep_CO_F", dir="forward", src=f"{EP}/path2c_CO-extrusion-TS_forward",
                    cc=[f"{CC2}/path2c_CO-extrusion-TS_forward/output.out"], classify="co",
                    label="CO-extrusion IRC end (forward)"),
               dict(id="ep_CO_B", dir="backward", src=f"{EP}/path2c_CO-extrusion-TS_backward",
                    cc=[f"{CC2}/path2c_CO-extrusion-TS_backward/output.out"], classify="co",
                    label="CO-extrusion IRC end (backward)")]),
]

# ------------------------------------------------------------------------------------------------ other paths
# Trajectories only (frames with their energies), shown in the Paths tab after the IRCs above.
PATHS = [
    dict(id="irc_a1_CO", label="IRC: CO extrusion (attempt 1)", group="Earlier IRCs (attempt 1)", keys=K_CO, ref="R",
         src=["attempt1/orca6_IRC/67-alkyne-tetracyclone-decarbonylation-TS/orca_IRC_Full_trj.xyz"],
         note="IRC from the CO-extrusion TS with the template's loose thresholds: one end is the [4+2] adduct, the other "
              "the decarbonylated product with CO 3.2 Å away; both ends stop before the bottom of their wells."),
    dict(id="irc_a1_DA", label="IRC: [4+2], higher TS (attempt 1)", group="Earlier IRCs (attempt 1)", keys=K_DA, ref="R",
         src=["attempt1/orca6_IRC/67-alkyne-tetracyclone-trap-TS/orca_IRC_Full_trj.xyz"],
         note="IRC from the attempt-1 [4+2] TS (the higher conformer). The adduct side is complete; the reactant side "
              "stops 3–4 kJ/mol below the TS on the flat entrance channel."),
    dict(id="irc_a1_SSE", label="IRC: S,S,E [1+2] (attempt 1)", group="Earlier IRCs (attempt 1)", keys=K_E, ref="R",
         src=["attempt1/orca6_IRC/67-carbene-tetracyclone-SSE-trap-TSb/orca_IRC_Full_trj.xyz"],
         note="Superseded by the complete S,S,E IRC above. The cyclopropane side is complete; the reactant side stops "
              "3–4 kJ/mol below the TS."),
    dict(id="neb_FBW", label="NEB-TS band: FBW", group="NEB-TS bands (attempt 3)", keys=K_FBW, ref="carbene",
         src=["attempt3/Step2_NEB-TS/path2a_FBW-TS/orca_MEP_trj.xyz"],
         note="Converged NEB band from the carbene to the alkyne (12 images). Its TS step reproduced the original FBW TS."),
    dict(id="neb_SSZ", label="NEB-TS band: S,S,Z [1+2]", group="NEB-TS bands (attempt 3)", keys=K_Z, ref="R",
         src=["attempt3/Step2_NEB-TS/path1a_SSZ-TSa/orca_MEP_trj.xyz"],
         note="Converged band from the Z complex to the S,S,Z product. The band crosses about 20 kJ/mol above the TS (the "
              "carbene cage is turned); ORCA's TS step from its top image relaxed onto S,S,Z-TSa."),
    dict(id="neb_SSE", label="NEB-TS band: S,S,E [1+2]", group="NEB-TS bands (attempt 3)", keys=K_E, ref="R",
         src=["attempt3/Step2_NEB-TS/path1b_SSE-TSb/orca_MEP_trj.xyz"],
         note="Converged band from the E complex to the S,S,E product. As for S,S,Z, the band crosses higher than the TS, "
              "and ORCA's TS step relaxed onto S,S,E-TSb."),
    dict(id="scan_DA", label="Relaxed scan: [4+2]", group="Relaxed scans (attempt 2)", keys=K_DA, ref="R",
         src="attempt2/Step2_scan/67-alkyne-tetracyclone-trap-TS/orca.[0-9][0-9][0-9].xyz",
         note="Both new C–C bonds stretched together from the adduct (0.05 Å steps). The maximum (2.48 Å) matches the "
              "attempt-1 [4+2] TS; the lowest [4+2] TS, found later from the long NEB, is 2.5 kJ/mol lower."),
    dict(id="scan_CO", label="Relaxed scan: CO extrusion", group="Relaxed scans (attempt 2)", keys=K_CO, ref="R",
         src="attempt2/Step2_scan/67-alkyne-tetracyclone-decarbonylation-TS/orca.[0-9][0-9][0-9].xyz",
         note="Both C–CO bonds stretched together from the adduct. The maximum (1.96 Å) is 0.4 kJ/mol below the TS."),
    dict(id="scan_SSZ", label="Relaxed scan: S,S,Z [1+2]", group="Relaxed scans (attempt 2)", keys=K_Z, ref="R",
         src="attempt2/Step2_scan/67-carbene-tetracyclone-SSZ-trap-TSa/orca.[0-9][0-9][0-9].xyz",
         note="C15–C22 shortened from 3.5 Å. The maximum (2.25 Å) matches the TS energy to 0.04 kJ/mol."),
    dict(id="scan_SSE", label="Relaxed scan: S,S,E [1+2]", group="Relaxed scans (attempt 2)", keys=K_E, ref="R",
         src="attempt2/Step2_scan/67-carbene-tetracyclone-SSE-trap-TSb/orca.[0-9][0-9][0-9].xyz",
         note="C15–C21 shortened from 3.5 Å. The maximum (2.25 Å) is 1.4 kJ/mol above the TS."),
]

# ------------------------------------------------------------------------------------------------ diagram states
# state id -> (label, members, short label). Each member is a list of candidate structures: the first one that has
# frequencies is used, otherwise the first one that has an energy. So an IRC end point takes over from the older
# structure as soon as its frequencies are done ("ep_CO_adduct" and "ep_CO_prodCO" are the CO-IRC ends, assigned
# from their geometry).
STATES = {
    "R": ("[2.2.2] carbene + tetracyclone", [["carbene"], ["tetracyclone"]], "carbene + tetracyclone"),
    "cZ": ("Carbene·tetracyclone complex (Z)", [["ep_SSZ_F", "complex_Z"]], "complex (Z)"),
    "tsZ": ("S,S,Z-#-TSa", [["TS_SSZ"]], "S,S,Z-TSa"),
    "pZ": ("S,S,Z-#", [["ep_SSZ_B", "SSZ"]], "S,S,Z-#"),
    "cE": ("Carbene·tetracyclone complex (E)", [["ep_SSE_B", "complex_E"]], "complex (E)"),
    "tsE": ("S,S,E-#-TSb", [["TS_SSE"]], "S,S,E-TSb"),
    "pE": ("S,S,E-#", [["ep_SSE_F", "SSE"]], "S,S,E-#"),
    "carb": ("[2.2.2] carbene", [["carbene"]], "carbene"),
    "tsFo": ("#-FBW-TS", [["TS_FBW"]], "FBW-TS"),
    "alko": ("[3.2.2] alkyne", [["alkyne"]], "alkyne"),
    "tsF": ("#-FBW-TS + tetracyclone", [["TS_FBW"], ["tetracyclone"]], "FBW-TS"),
    "alk": ("[3.2.2] alkyne + tetracyclone", [["alkyne"], ["tetracyclone"]], "alkyne + tetracyclone"),
    "cA": ("Alkyne·tetracyclone complex", [["ep_DA_B", "complex_alkyne"]], "complex"),
    "tsDA": ("#-TS [4+2]", [["TS_DA"]], "TS [4+2]"),
    "tsDA1": ("#-TS [4+2], higher conformer", [["TS_DA_a1"]], "TS [4+2] (a1)"),
    "addB": ("[4+2] adduct from the lowest [4+2] TS", [["ep_DA_F"]], "adduct (IRC)"),
    "add": ("[4+2] adduct", [["ep_CO_adduct", "adduct"]], "[4+2] adduct"),
    "tsCO": ("#-TS CO extrusion", [["TS_CO"]], "TS CO extrusion"),
    "cP": ("Decarbonylated product·CO complex", [["ep_CO_prodCO", "complex_prodCO"]], "product·CO complex"),
    "prod": ("Decarbonylated product + CO", [["product"], ["CO"]], "product + CO"),
}

# Routes drawn as energy diagrams. Links: (from, to) or (from, to, "conf") for a conformational change that has not
# been computed (drawn dotted).
ROUTES = [
    dict(id="overview", label="Whole surface (as in the PDF)",
         desc="Section B as drawn in the PDF: the two [1+2] trappings on the left, the FBW rearrangement, Diels–Alder "
              "and CO extrusion on the right. Pre-reaction complexes are left out here, as in the PDF.",
         states=["pZ", "pE", "tsZ", "tsE", "R", "tsF", "alk", "tsDA", "add", "tsCO", "prod"],
         links=[("pZ", "tsZ"), ("tsZ", "R"), ("pE", "tsE"), ("tsE", "R"), ("R", "tsF"), ("tsF", "alk"),
                ("alk", "tsDA"), ("tsDA", "add"), ("add", "tsCO"), ("tsCO", "prod")]),
    dict(id="route1a", label="Path 1a: [1+2] trapping, Z face",
         desc="The carbene adds across a C=C bond of tetracyclone on the Z face to give the S,S,Z cyclopropane.",
         states=["R", "cZ", "tsZ", "pZ"], links=[("R", "cZ"), ("cZ", "tsZ"), ("tsZ", "pZ")]),
    dict(id="route1b", label="Path 1b: [1+2] trapping, E face",
         desc="The same addition on the E face, giving the S,S,E cyclopropane.",
         states=["R", "cE", "tsE", "pE"], links=[("R", "cE"), ("cE", "tsE"), ("tsE", "pE")]),
    dict(id="route2", label="Path 2: FBW → [4+2] → CO extrusion",
         desc="The carbene first rearranges to the strained [3.2.2] alkyne (FBW), which tetracyclone traps in a "
              "Diels–Alder reaction; the adduct then loses CO. The lowest [4+2] TS leads to a different conformer of the "
              "adduct from the one the CO-extrusion TS starts from (dotted link: their interconversion has not been "
              "computed).",
         states=["R", "tsF", "alk", "cA", "tsDA", "addB", "add", "tsCO", "cP", "prod"],
         links=[("R", "tsF"), ("tsF", "alk"), ("alk", "cA"), ("cA", "tsDA"), ("tsDA", "addB"), ("addB", "add", "conf"),
                ("add", "tsCO"), ("tsCO", "cP"), ("cP", "prod")]),
]

# Relative energies grouped by reaction (as in the SI): every row relative to the first (the separated reactants).
REACTIONS = [
    dict(id="fbw", label="FBW rearrangement, 6→7 ([2.2.2] carbene → [3.2.2] alkyne)", rows=["carb", "tsFo", "alko"]),
    dict(id="z12", label="[1+2] trapping by tetracyclone, Z face", rows=["R", "cZ", "tsZ", "pZ"]),
    dict(id="e12", label="[1+2] trapping by tetracyclone, E face", rows=["R", "cE", "tsE", "pE"]),
    dict(id="da", label="[4+2] trapping of the [3.2.2] alkyne by tetracyclone",
         rows=["alk", "cA", "tsDA", "addB", "tsDA1", "add"]),
    dict(id="co", label="CO extrusion from the [4+2] adduct", rows=["add", "tsCO", "cP", "prod"]),
    dict(id="all", label="Whole surface, relative to [2.2.2] carbene + tetracyclone (the PDF's 0)",
         rows=["R", "cZ", "tsZ", "pZ", "cE", "tsE", "pE", "tsF", "alk", "cA", "tsDA", "addB", "add", "tsCO", "cP",
               "prod"]),
]
