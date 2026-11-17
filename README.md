# FBW tetracyclone trapping explorer

An interactive viewer for the reactants, intermediates, transition states and products of the 6→7
Fritsch–Buttenberg–Wiechell (FBW) system trapped by tetracyclone (section B of the study's ORCA work summary).
The [2.2.2] carbene either adds directly to tetracyclone ([1+2], two faces), or first rearranges to the [3.2.2] alkyne
(FBW), which tetracyclone traps by a Diels–Alder reaction followed by CO extrusion.

**Live site:** https://geoffreyweal.github.io/FBW_Tetracyclone_Trapping_Explorer/. You can also open `index.html` directly in a browser: no server or installation is
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

- **Structures:** M06-2X/def2-TZVP D3(0), CPCM(CH₂Cl₂: ε = 9.08, n = 1.424), RIJCOSX (def2/J), DefGrid2, TightSCF, TightOPT; ORCA 6.1.0. Minima have no imaginary frequency and transition structures exactly one.
  `thermochemistry/method_check.csv` compares every input with this.
- **Thermochemistry:** GoodVibes 4.3.0 at 313.15 K from the M06-2X harmonic frequencies: Grimme quasi-RRHO entropy (100 cm⁻¹ cut-off), no quasi-harmonic enthalpy correction, no frequency scaling; 1 mol/L standard state. The 1 mol/L values add RT ln(Vm/V°) = 8.45 kJ/mol to the
  1 atm value of each species. Commands (run by `tools/build_site.py`):
  - `goodvibes --temp 313.15 -v 1.0 --imag --check -c 1.0 *.out`
  - `goodvibes --temp 313.15 -v 1.0 --imag --check *.out` (1 atm values)
- **Single points:** DLPNO-CCSD(T)/def2-TZVPP (def2-TZVPP/C, TightPNO), same CPCM, on the M06-2X geometries. Composite G = Eel(DLPNO-CCSD(T)) + [G − Eel](M06-2X), from GoodVibes `--spc`.
  4 of the 19 structures in the diagrams have their DLPNO-CCSD(T) energy so far, so
  **all relative energies below are M06-2X values** for now.
- **Units:** kJ/mol, 1 Eh = 2625.50 kJ/mol, relative to the separated reactants of each step.

## Relative energies by reaction (kJ/mol, M06-2X, 313.15 K)

**FBW rearrangement, 6→7 ([2.2.2] carbene → [3.2.2] alkyne)**

| Structure | ΔEel | ΔG (1 M) | ΔG (1 atm) |
| --- | --- | --- | --- |
| [2.2.2] carbene | +0.0 | +0.0 | +0.0 |
| #-FBW-TS | +28.1 | +29.7 | +29.7 |
| [3.2.2] alkyne | -50.3 | -44.7 | -44.7 |

**[1+2] trapping by tetracyclone, Z face**

| Structure | ΔEel | ΔG (1 M) | ΔG (1 atm) |
| --- | --- | --- | --- |
| [2.2.2] carbene + tetracyclone | +0.0 | +0.0 | +0.0 |
| Carbene·tetracyclone complex (Z) | -38.8 | +13.0 | +21.4 |
| S,S,Z-#-TSa | -17.8 | +45.4 | +53.9 |
| S,S,Z-# | -322.1 | -250.6 | -242.1 |

**[1+2] trapping by tetracyclone, E face**

| Structure | ΔEel | ΔG (1 M) | ΔG (1 atm) |
| --- | --- | --- | --- |
| [2.2.2] carbene + tetracyclone | +0.0 | +0.0 | +0.0 |
| Carbene·tetracyclone complex (E) | -45.4 | +10.8 | +19.2 |
| S,S,E-#-TSb | -27.1 | +33.8 | +42.2 |
| S,S,E-# | -322.4 | -250.0 | -241.5 |

**[4+2] trapping of the [3.2.2] alkyne by tetracyclone**

| Structure | ΔEel | ΔG (1 M) | ΔG (1 atm) |
| --- | --- | --- | --- |
| [3.2.2] alkyne + tetracyclone | +0.0 | +0.0 | +0.0 |
| Alkyne·tetracyclone complex | -39.6 | – | – |
| #-TS [4+2] | -16.2 | +42.1 | +50.5 |
| [4+2] adduct from the lowest [4+2] TS | -278.4 | – | – |
| #-TS [4+2], higher conformer | -13.7 | +43.9 | +52.4 |
| [4+2] adduct | -279.1 | -206.6 | -198.1 |

**CO extrusion from the [4+2] adduct**

| Structure | ΔEel | ΔG (1 M) | ΔG (1 atm) |
| --- | --- | --- | --- |
| [4+2] adduct | +0.0 | +0.0 | +0.0 |
| #-TS CO extrusion | +79.5 | +71.9 | +71.9 |
| Decarbonylated product·CO complex | -134.7 | – | – |
| Decarbonylated product + CO | -114.1 | -166.3 | -174.7 |

## Whole surface, relative to [2.2.2] carbene + tetracyclone (the PDF's 0)

| State | ΔEel | ΔG (1 M) | ΔG (1 atm) | Structures used |
| --- | --- | --- | --- | --- |
| [2.2.2] carbene + tetracyclone | +0.0 | +0.0 | +0.0 | `carbene` + `tetracyclone` |
| Carbene·tetracyclone complex (Z) | -38.8 | +13.0 | +21.4 | `complex_Z` |
| S,S,Z-#-TSa | -17.8 | +45.4 | +53.9 | `TS_SSZ` |
| S,S,Z-# | -322.1 | -250.6 | -242.1 | `ep_SSZ_B` |
| Carbene·tetracyclone complex (E) | -45.4 | +10.8 | +19.2 | `ep_SSE_B` |
| S,S,E-#-TSb | -27.1 | +33.8 | +42.2 | `TS_SSE` |
| S,S,E-# | -322.4 | -250.0 | -241.5 | `ep_SSE_F` |
| #-FBW-TS + tetracyclone | +28.1 | +29.7 | +29.7 | `TS_FBW` + `tetracyclone` |
| [3.2.2] alkyne + tetracyclone | -50.3 | -44.7 | -44.7 | `alkyne` + `tetracyclone` |
| Alkyne·tetracyclone complex | -89.8 | – | – | `ep_DA_B` |
| #-TS [4+2] | -66.5 | -2.7 | +5.8 | `TS_DA` |
| [4+2] adduct from the lowest [4+2] TS | -328.7 | – | – | `ep_DA_F` |
| [4+2] adduct | -329.4 | -251.3 | -242.8 | `adduct` |
| #-TS CO extrusion | -250.0 | -179.4 | -171.0 | `TS_CO` |
| Decarbonylated product·CO complex | -464.1 | – | – | `complex_prodCO` |
| Decarbonylated product + CO | -443.5 | -417.6 | -417.6 | `product` + `CO` |

– means no frequencies yet. An IRC end point is used as soon as its frequencies are done; until then the earlier
structure (or the end point's electronic energy alone) stands in.

## Transition structures

| Transition structure | ν_imag (cm⁻¹) | Forming (Å) | Breaking (Å) |
| --- | --- | --- | --- |
| S,S,Z-#-TSa ([1+2], Z face) | -214.0 | 2.242, 2.611 | — |
| S,S,E-#-TSb ([1+2], E face) | -230.3 | 2.221, 2.600 | — |
| #-FBW-TS (FBW rearrangement) | -362.5 | 1.797 | 1.970 |
| #-TS [4+2] (Diels–Alder) | -237.6 | 2.461, 2.465 | — |
| #-TS [4+2], higher conformer (attempt 1) | -212.5 | 2.496, 2.471 | — |
| #-TS CO extrusion | -446.6 | — | 1.987, 1.974 |

FBW shift: distance from the migrating carbon (C5) to its destination (C16, forming) and origin (C2, breaking).
Additions: the two C···C distances between the reacting fragments. CO extrusion: the two breaking C–C(O) bonds.

## IRCs and their optimised end points

| IRC | Direction | End point | State | Nimag | Compared with |
| --- | --- | --- | --- | --- | --- |
| IRC: S,S,Z [1+2] | forward | Carbene·tetracyclone complex, Z face (S,S,Z IRC end) | freq |  | same as `complex_Z` (RMSD 0.233 Å, ΔE +0.10 kJ/mol) |
| IRC: S,S,Z [1+2] | backward | S,S,Z-# (S,S,Z IRC end) | complete | 0 | same as `SSZ` (RMSD 0.127 Å, ΔE -0.02 kJ/mol, ΔG +1.24) |
| IRC: S,S,E [1+2] | forward | S,S,E-# (S,S,E IRC end) | complete | 0 | same as `SSE` (RMSD 0.010 Å, ΔE -0.04 kJ/mol, ΔG +0.24) |
| IRC: S,S,E [1+2] | backward | Carbene·tetracyclone complex, E face (S,S,E IRC end) | complete | 0 | same as `complex_E` (RMSD 0.083 Å, ΔE +0.09 kJ/mol, ΔG +0.48) |
| IRC: FBW rearrangement | forward | [2.2.2] carbene (FBW IRC end) | complete | 0 | same as `carbene` (RMSD 0.001 Å, ΔE -0.00 kJ/mol, ΔG +0.01) |
| IRC: FBW rearrangement | backward | [3.2.2] alkyne (FBW IRC end) | complete | 0 | same as `alkyne` (distance RMS (atom order differs) 0.013 Å, ΔE -0.17 kJ/mol, ΔG -0.41) |
| IRC: [4+2] (lowest TS) | forward | [4+2] adduct from the lowest [4+2] TS ([4+2] IRC end) | freq |  | differs from `adduct` (RMSD 0.677 Å, ΔE +0.71 kJ/mol) |
| IRC: [4+2] (lowest TS) | backward | Alkyne·tetracyclone complex ([4+2] IRC end) | freq |  | same as `complex_alkyne` (RMSD 0.003 Å, ΔE -0.00 kJ/mol) |
| IRC: CO extrusion | forward | – | not started |  | – |
| IRC: CO extrusion | backward | – | not started |  | – |

## Still to come

- Carbene·tetracyclone complex, Z face (S,S,Z IRC end): optimised; frequencies still running
- [4+2] adduct from the lowest [4+2] TS ([4+2] IRC end): optimised; frequencies still running
- Alkyne·tetracyclone complex ([4+2] IRC end): optimised; frequencies still running
- IRC: CO extrusion: the IRC is still running; Opt + Freq of its two ends follow
- No free energy yet (electronic energy only) for: Alkyne·tetracyclone complex; [4+2] adduct from the lowest [4+2] TS; Decarbonylated product·CO complex
- DLPNO-CCSD(T) single points: 4 of the 19 structures in the diagrams are done; relative energies stay at the M06-2X level until all are

## Absolute energies (Eh)

| Structure | Type | Nimag | Eel M06-2X | G − Eel | G (M06-2X, 1 M) | Eel DLPNO-CCSD(T) | T1 | File |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| [2.2.2] carbene | minimum | 0 | -350.022495 | 0.156505 | -349.865990 | -349.410160 | 0.012 | `xyz/species/carbene.xyz` |
| Tetracyclone (trapping reagent) | minimum | 0 | -1192.295058 | 0.351403 | -1191.943655 | -1190.165957 | 0.011 | `xyz/species/tetracyclone.xyz` |
| Carbene·tetracyclone complex, Z face | complex | 0 | -1542.332343 | 0.527642 | -1541.804701 | – |  | `xyz/species/complex_Z.xyz` |
| S,S,Z-#-TSa ([1+2], Z face) | ts | 1 | -1542.324348 | 0.532009 | -1541.792339 | – |  | `xyz/species/TS_SSZ.xyz` |
| S,S,Z-# ([1+2] adduct) | minimum | 0 | -1542.440208 | 0.534664 | -1541.905544 | – |  | `xyz/species/SSZ.xyz` |
| Carbene·tetracyclone complex, E face | complex | 0 | -1542.334880 | 0.529167 | -1541.805714 | – |  | `xyz/species/complex_E.xyz` |
| S,S,E-#-TSb ([1+2], E face) | ts | 1 | -1542.327858 | 0.531072 | -1541.796786 | – |  | `xyz/species/TS_SSE.xyz` |
| S,S,E-# ([1+2] adduct) | minimum | 0 | -1542.440333 | 0.535388 | -1541.904945 | – |  | `xyz/species/SSE.xyz` |
| #-FBW-TS (FBW rearrangement) | ts | 1 | -350.011791 | 0.157109 | -349.854682 | -349.395114 | 0.011 | `xyz/species/TS_FBW.xyz` |
| [3.2.2] alkyne | minimum | 0 | -350.041641 | 0.158612 | -349.883029 | -349.429325 | 0.011 | `xyz/species/alkyne.xyz` |
| Alkyne·tetracyclone complex (attempt 3) | complex |  | -1542.351766 | – | – | – |  | `xyz/species/complex_alkyne.xyz` |
| #-TS [4+2] (Diels–Alder) | ts | 1 | -1542.342879 | 0.532213 | -1541.810666 | – |  | `xyz/species/TS_DA.xyz` |
| #-TS [4+2], higher conformer (attempt 1) | ts | 1 | -1542.341926 | 0.531961 | -1541.809965 | – |  | `xyz/species/TS_DA_a1.xyz` |
| [4+2] adduct | minimum | 0 | -1542.443021 | 0.537663 | -1541.905358 | – |  | `xyz/species/adduct.xyz` |
| #-TS CO extrusion | ts | 1 | -1542.412761 | 0.534780 | -1541.877981 | – |  | `xyz/species/TS_CO.xyz` |
| Decarbonylated product·CO complex (attempt 3) | complex |  | -1542.494313 | – | – | – |  | `xyz/species/complex_prodCO.xyz` |
| Decarbonylated product | minimum | 0 | -1429.163521 | 0.529643 | -1428.633878 | – |  | `xyz/species/product.xyz` |
| Carbon monoxide | minimum | 0 | -113.322954 | -0.011853 | -113.334808 | – |  | `xyz/species/CO.xyz` |
| Carbene·tetracyclone complex, Z face (S,S,Z IRC end) | complex |  | -1542.332303 | – | – | – |  | `xyz/species/ep_SSZ_F.xyz` |
| S,S,Z-# (S,S,Z IRC end) | minimum | 0 | -1542.440214 | 0.535141 | -1541.905074 | – |  | `xyz/species/ep_SSZ_B.xyz` |
| S,S,E-# (S,S,E IRC end) | minimum | 0 | -1542.440347 | 0.535494 | -1541.904853 | – |  | `xyz/species/ep_SSE_F.xyz` |
| Carbene·tetracyclone complex, E face (S,S,E IRC end) | complex | 0 | -1542.334848 | 0.529317 | -1541.805531 | – |  | `xyz/species/ep_SSE_B.xyz` |
| [2.2.2] carbene (FBW IRC end) | minimum | 0 | -350.022495 | 0.156509 | -349.865987 | – |  | `xyz/species/ep_FBW_F.xyz` |
| [3.2.2] alkyne (FBW IRC end) | minimum | 0 | -350.041706 | 0.158520 | -349.883186 | – |  | `xyz/species/ep_FBW_B.xyz` |
| [4+2] adduct from the lowest [4+2] TS ([4+2] IRC end) | minimum |  | -1542.442750 | – | – | – |  | `xyz/species/ep_DA_F.xyz` |
| Alkyne·tetracyclone complex ([4+2] IRC end) | complex |  | -1542.351767 | – | – | – |  | `xyz/species/ep_DA_B.xyz` |

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

Built 2026-11-18.
