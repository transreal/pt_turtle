# pt_turtle — turtle walks on Penrose rhomb tilings

Computer-assisted proofs, machine-checkable certificates and exhaustive surveys for the
boundedness of *turtle* walks on Penrose rhomb tilings.  This repository is the supplementary
material of the paper

> Katsunobu Imai, *Boundedness of turtle orbits on Penrose rhomb tilings* (prepared for ASCAT 2027).

The paper itself is not included before publication; the proofs are given here in
`proof_sketch.md` and `proof_unbounded.md`.  The question comes from the
reversible self-reproducing cellular automaton SR8 (Morita and Imai, 1996) carried over to Penrose
tilings in the preceding paper *Signal propagation and reversibility on Penrose tiling* (in:
Understanding Reversibility, World Scientific, to appear): the head of an SR8 worm moves like a
turtle, and the shortest worm follows the program (SSR)*.

## Turtles

A state is a pair (tile, edge) of a Penrose rhomb tiling; the edge is the heading.  The commands are
`S` (cross the heading edge into the neighbouring rhomb and take the opposite edge as the new heading,
which continues straight along a de Bruijn ribbon), `R` and `L` (rotate the heading clockwise /
counter-clockwise inside the tile).  A program P is applied cyclically, written (P)*.  Since S, R
and L are bijections of the set of states, every orbit is either periodic or infinite; a program is
*bounded* if every orbit on every Penrose rhomb tiling is periodic.

## Main results

* **(SSR)\* is bounded on every Penrose rhomb tiling.**  The periods are exactly
  {10, 12, 28, 30, 50, 80, 100} and the swept sets fall into 7 congruence classes.  Hence the
  shortest SR8 worm never collides with itself wherever it is placed on a Penrose tiling (its
  periods are 8 times the turtle periods).
* **59 programs are bounded** (52 up to mirror, rotation and powers), e.g. (SSR)(SSSSR),
  (SR)^3(SSR), (SSR)(SR)^3(SSR), and the single-exit family (SR)^{k-3}(SSR) for
  k = 3, 6, 9, 10, 12, 15, 18, 20, 21, 22, 24, 27, 29, 30, 35, 40, 45, 50.
  Method: exact Robinson-triangle substitution σ in Z[ζ], an atlas of 30 marked-triangle coronas
  that is closed under σ (hence complete, by repetitivity), and exhaustive verification of all
  orbits starting in the central supertriangle of the patches σ^n(corona), n = 7 … 10, none of
  which attempts to leave its patch.  The verification was repeated with three independent
  implementations (the Python code here, `wl_verify.wl` on top of the author's `TurtleTiling.wl`,
  and the Wolfram Language code in `notebook_src/`).
* **(SSSR)\* is not bounded.**  In every Penrose rhomb tiling and for every odd m ≥ 7, every
  positively oriented obtuse level-m supertriangle carries an orbit segment of
  N_m = 7·2^(m−2) − 22 S-moves (202, 874, 3562, …) with pairwise distinct states; so periods are
  unbounded in every tiling, and by compactness some Penrose tiling has an infinite orbit.  The
  proof is an induction on the substitution level over a finite *pass-substitution system*
  (60 pass families and 120 port families obeying the exact affine law a_{m+2} = φ² a_m + δ),
  given as a certificate and accepted by two independently written checkers.  The same method
  proves 18 programs unbounded, e.g. (SLSR)*, (SR)(SSR), (SR)^2(SSR), (SSSSSR)*.
* **Three regimes.**  On windows of side 46 / 160 / 300 (pentagrid coordinates) of the pentagrid
  tiling with offsets (0, 0.1, −0.1, 0.2, −0.2), the 43 programs of length ≤ 8 without adjacent
  turns split into uniformly bounded (7), hierarchical (30; orbit sizes grow with the window) and
  diffusive (6; e.g. (SSLSR)*, almost all orbits escape, displacement ~ t^0.4 … t^0.5).

## Layout

| path | content |
|---|---|
| `proof_sketch.md` | the boundedness argument for (SSR)* (English) |
| `proof_unbounded.md` | statement, proof and certificate format of the unboundedness of (SSSR)* (English) |
| `pt.py` | pentagrid tiling (de Bruijn), edge adjacency, vectorised turtle, orbit classification |
| `subst.py`, `subst_fast.py` | exact Robinson-triangle substitution in Z[ζ], coronas, canonical forms |
| `prove2.py`, `prove_many.py` | boundedness pipeline: atlas, closure check, verification V(P, n), certificates |
| `export_patches.py`, `wl_verify.wl` | export of the verification patches and an independent Wolfram Language check (needs the author's `CellularAutomata.wl` / `TurtleTiling.wl`, not included) |
| `turtle_survey.cu`, `build.bat`, `survey.py`, `survey_all.py`, `survey_exits.py`, `survey_ronly.py`, `survey_ronly2.py`, `progs.py` | CUDA survey kernel and drivers of the window surveys |
| `test_pivot.py`, `diffusion.py`, `figs.py`, `fig_pass.py`, `tables.py` | pivot-walk check, diffusion exponent, figures (`figs/`), tables |
| `lsys.py`, `extract2.py`, `batch_unbounded.py`, `collect_certs.py` | extraction of the pass-substitution certificates |
| `verify_cert.py`, `indep_check.py` | the two independent certificate checkers |
| `results/` | boundedness certificates (`proof2_*.json`, `proof_many.json`), survey outputs (`surveyall_*.json`, `surveyexits_300.json`, `surveyronly*.json`), unboundedness certificates (`unbounded_certs/`) |
| `logs/`, `figs/` | logs of the runs and figures produced by the Python scripts |
| `notebook_src/` | Wolfram Language code of a third, self-contained implementation (exact substitution, atlas, verification, exact congruence classes, window surveys) and the figure scripts of the paper; see `notebook_src/README.md` |

## Reproduce

```
python prove2.py SSR 7 8                 # (SSR)*: atlas, closure and verification at level 7 (about 2 min)
python prove2.py SSRSSSSR 9 8 --cuda     # (SSR)(SSSSR) at level 9 on the GPU
python prove_many.py --auto              # all bounded programs of results/proof_many.json
python survey.py gen 300                 # window of side 300, then e.g.
python survey.py run 300 SSSR            # survey of (SSSR)* (GPU kernel built with build.bat)
python verify_cert.py results/unbounded_certs/cert2_SSSR_1p_1.json     # checker 1
python indep_check.py results/unbounded_certs/cert2_SSSR_1p_1.json     # checker 2
python export_patches.py SSR 7           # writes results/patches_SSR_n7.json for wl_verify.wl
wolframscript -file wl_verify.wl results/patches_SSR_n7.json
```

Python 3.14 with numpy 2.4 (matplotlib for the figures); CUDA 12.9 and MSVC for `turtle_survey.cu`;
Wolfram 15.0 for `notebook_src/` and `wl_verify.wl`.

## Notes

* The orbit-class counts written by the GPU survey (`results/surveyall_*.json`, `surveyexits_300.json`,
  `surveyronly*.json`) come from a float32 signature and are too large for programs with large orbits
  (e.g. (SSRSSSSR)*: 38 instead of 13; (SSSR)* on the window of side 300: 63 instead of 41).  The paper
  uses exact congruence classes computed in `notebook_src/` (`survey_exact.wl`); the counts of closed
  and escaped orbits agree exactly.
* `results/patches_SSR_n7.json` (24 MB, the input of `wl_verify.wl`) is not included; `export_patches.py`
  regenerates it.

## Citation

Please cite the ASCAT 2027 paper, and the preceding paper for SR8 on Penrose tilings:
Katsunobu Imai, *Signal propagation and reversibility on Penrose tiling*, in: Understanding
Reversibility, World Scientific (to appear).

## Revision of 2026-10-07 (tag `ascat2027`)

Changes made after a critical review of the ASCAT 2027 manuscript:

* `PROGRAMS.md` — complete lists of the programs of the theorems: the 59 bounded programs (with the
  verification level, the number of initial states, the period spectrum, the number of swept-set classes
  and the agreement of the two implementations), the 15 further bounded programs, and the 21
  unboundedness certificates with their root pass family, lengths, exact recurrence and checker results.
* `notebook_src/verify_orbit_classes.wl` — exact classification of the closed orbits of a program on the
  verification patches as *cyclic sequences of states* (rhomb vertices and heading edge in Z[zeta]) up to
  translations and rotations by multiples of 36 degrees, i.e. the orientation-preserving part of G.
  For (SSR)* at level 7 this gives 7 orbits, one for each period, so the 7 swept-set classes of the paper
  are also 7 dynamical orbit classes (log in `notebook_src/data/orbit_classes_SSR_n7.log`).
* `notebook_src/lib.wl`, `verify_classes.wl` — congruence classes are identified by comparing the exact
  canonical forms themselves; the earlier version compared hash values of the canonical forms.
* `notebook_src/review_data.wl` — the numbers quoted in the survey section: initial states per window
  (tiles whose four neighbours lie in the window: 15 902 / 195 801 / 690 036 tiles, i.e. 63 608 / 783 204 /
  2 760 144 states for the windows 46 / 160 / 300) and the distance from the centre of each of the 30 coronas
  to the boundary of the corona (minimum sin 36 degrees, up to rounding).
* Exactness: all coincidences of points and all congruence tests are integer computations in Z[zeta].
  The only floating-point quantities in the proofs are the distance margins of obligation O5 / (V3) of the
  pass-substitution certificates, whose minimum margin (1.89 for (SSSR)*) exceeds the rounding error of the
  double-precision evaluation (below 1e-9) by many orders of magnitude.
* The survey of the "43 programs" enumerates the programs of length at most 8 that contain `S` and at least
  one turn, with no two turns cyclically adjacent, up to cyclic shifts and mirroring (powers are not
  identified); `progs.py`.
