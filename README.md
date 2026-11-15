# FBW tetracyclone trapping explorer

An interactive viewer for the reactants, intermediates, transition states and products of the 6→7
Fritsch–Buttenberg–Wiechell (FBW) system trapped by tetracyclone. This is section B of the study's ORCA work summary:
the [2.2.2] carbene either adds directly to tetracyclone ([1+2], two faces), or first rearranges to the [3.2.2] alkyne
(FBW), which tetracyclone traps by a Diels–Alder reaction followed by CO extrusion.

**Live site:** https://geoffreyweal.github.io/FBW_Tetracyclone_Trapping_Explorer/ (once GitHub Pages is switched on).
You can also open `index.html` directly in a browser: no server or installation is needed, and the page makes no
network connections.

## What you can do

- **Routes:** energy diagrams of each route (qh-G or electronic energy). Click a level to see its structure in 3D, with
  the bonds that form or break shown as dashed lines with their lengths.
- **Species:** every structure in a table with its energies and imaginary frequency; click a row to view it.
- **Paths:** step through IRCs, NEB bands and relaxed scans frame by frame, with the energy along the path.
- **Download:** every structure is in `xyz/` (open with `ase gui`, Avogadro, Jmol and so on).

## Energies

M06-2X/def2-TZVP D3Zero, CPCM (ε = 9.08, n = 1.424), defgrid2, TightSCF, RIJCOSX; ORCA 6.1.0. Quasi-RRHO free energies (qh-G) as in GoodVibes: Grimme entropy with a 100 cm⁻¹ cut-off, 313.15 K, 1 mol/L standard state. All energies in kJ/mol relative to the [2.2.2] carbene and tetracyclone as separate molecules (the PDF's 0).

| State | E (kJ/mol) | qh-G (kJ/mol) |
| --- | --- | --- |
| [2.2.2] carbene + tetracyclone | +0.0 | +0.0 |
| Carbene·tetracyclone complex (Z) | -38.8 | +13.0 |
| S,S,Z-#-TSa | -17.8 | +45.4 |
| S,S,Z-# | -322.0 | -251.8 |
| Carbene·tetracyclone complex (E) | -45.5 | +10.3 |
| S,S,E-#-TSb | -27.1 | +33.8 |
| S,S,E-# | -322.4 | -250.2 |
| #-FBW-TS + tetracyclone | +28.1 | +29.7 |
| [3.2.2] alkyne + tetracyclone | -50.3 | -44.7 |
| Alkyne·tetracyclone complex | -89.8 | – |
| #-TS [4+2] | -64.0 | -0.8 |
| [4+2] adduct | -329.4 | -251.3 |
| #-TS CO extrusion | -250.0 | -179.4 |
| Decarbonylated product·CO complex | -464.1 | – |
| Decarbonylated product + CO | -443.5 | -417.6 |

– means no frequencies are available for that structure.

## Structures

| Structure | Type | Atoms | Imaginary frequency (cm⁻¹) | File | Source |
| --- | --- | --- | --- | --- | --- |
| [2.2.2] carbene | minimum | 21 | – | `xyz/species/carbene.xyz` | Original study (DFT minimum) |
| Tetracyclone (trapping reagent) | minimum | 50 | – | `xyz/species/tetracyclone.xyz` | Original study (DFT minimum; last step of its optimisation) |
| Carbene·tetracyclone complex, Z face | complex | 71 | – | `xyz/species/complex_Z.xyz` | Attempt 3, Step 1 (Opt + Freq, no imaginary frequencies) |
| S,S,Z-#-TSa ([1+2], Z face) | ts | 71 | -214.0 | `xyz/species/TS_SSZ.xyz` | Attempt 1 (OptTS, analytic Hessian, v2); reproduced by attempt-3 NEB-TS and the attempt-2 scan |
| S,S,Z-# ([1+2] adduct) | minimum | 71 | – | `xyz/species/SSZ.xyz` | Original study (DFT minimum) |
| Carbene·tetracyclone complex, E face | complex | 71 | – | `xyz/species/complex_E.xyz` | Attempt 3, Step 1 (Opt + Freq, no imaginary frequencies) |
| S,S,E-#-TSb ([1+2], E face) | ts | 71 | -230.3 | `xyz/species/TS_SSE.xyz` | Attempt 1 (OptTS, analytic Hessian, v2); reproduced by attempt-3 NEB-TS and the attempt-2 scan |
| S,S,E-# ([1+2] adduct) | minimum | 71 | – | `xyz/species/SSE.xyz` | Original study (DFT minimum) |
| #-FBW-TS (FBW rearrangement) | ts | 21 | -362.5 | `xyz/species/TS_FBW.xyz` | Original study; reproduced by attempt-3 NEB-TS (0.0002 Å, 0.02 kJ/mol) |
| [3.2.2] alkyne | minimum | 21 | – | `xyz/species/alkyne.xyz` | Original study (DFT minimum) |
| Alkyne·tetracyclone complex | complex | 71 | – | `xyz/species/complex_alkyne.xyz` | Attempt 3, Step 1 (optimisation converged; stopped before the frequencies) |
| #-TS [4+2] (Diels–Alder) | ts | 71 | -212.5 | `xyz/species/TS_DA.xyz` | Attempt 1 (OptTS, analytic Hessian, v2); matched by the attempt-2 scan maximum |
| [4+2] adduct | minimum | 71 | – | `xyz/species/adduct.xyz` | Original study (DFT minimum) |
| #-TS CO extrusion | ts | 71 | -446.6 | `xyz/species/TS_CO.xyz` | Attempt 1 (OptTS, analytic Hessian, v2); matched by the attempt-2 scan maximum |
| Decarbonylated product·CO complex | complex | 71 | – | `xyz/species/complex_prodCO.xyz` | Attempt 3, Step 1 (stopped at cycle 53, nearly converged: max gradient 2.6e-4; no frequencies) |
| Decarbonylated product | minimum | 69 | – | `xyz/species/product.xyz` | Original study (DFT minimum) |
| Carbon monoxide | minimum | 2 | – | `xyz/species/CO.xyz` | Original study (DFT minimum) |

## Paths

| Path | From | Frames | File |
| --- | --- | --- | --- |
| IRC: CO extrusion | IRC (attempt 1) | 116 | `xyz/paths/irc_CO.xyz` |
| IRC: [4+2] | IRC (attempt 1) | 48 | `xyz/paths/irc_DA.xyz` |
| IRC: S,S,E [1+2] | IRC (attempt 1) | 43 | `xyz/paths/irc_SSE.xyz` |
| NEB-TS band: FBW | NEB-TS bands (attempt 3) | 14 | `xyz/paths/neb_FBW.xyz` |
| NEB-TS band: S,S,Z [1+2] | NEB-TS bands (attempt 3) | 14 | `xyz/paths/neb_SSZ.xyz` |
| NEB-TS band: S,S,E [1+2] | NEB-TS bands (attempt 3) | 14 | `xyz/paths/neb_SSE.xyz` |
| Relaxed scan: [4+2] | Relaxed scans (attempt 2) | 32 | `xyz/paths/scan_DA.xyz` |
| Relaxed scan: CO extrusion | Relaxed scans (attempt 2) | 15 | `xyz/paths/scan_CO.xyz` |
| Relaxed scan: S,S,Z [1+2] | Relaxed scans (attempt 2) | 35 | `xyz/paths/scan_SSZ.xyz` |
| Relaxed scan: S,S,E [1+2] | Relaxed scans (attempt 2) | 37 | `xyz/paths/scan_SSE.xyz` |

Fresh IRCs from all five TSs, with both ends optimised to minima, are running and will be added when they finish.

## Files

- `index.html`, `assets/`: the viewer (3Dmol.js, BSD-3-Clause; see `assets/3Dmol_LICENSE.txt`).
- `data/`: the structures, energies and paths as small script files.
- `xyz/species/`, `xyz/paths/`: the same structures and paths as xyz files.

Built 2026-11-16 by `make_website.py` from the study's calculation outputs.
