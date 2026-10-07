# notebook_src — Wolfram Language verification code and figure scripts

A third, independent implementation of the computations of the Python code in the repository root
(exact Robinson-triangle substitution in Z[zeta] on the basis 1, zeta, zeta^2, zeta^3, corona atlas,
turtle, exact congruence classes of swept sets), together with the scripts that draw the figures of
the paper.  Everything here is Wolfram Language (15.0) except `py_extra.py`.  The text of the paper
is not included before publication.

| file | content |
|---|---|
| `lib.wl` | exact arithmetic, substitution, rhombs, navigation, compiled turtle, coronas, exact D10 canonical forms, pentagrid |
| `check_atlas.wl`, `verify_all.wl`, `verify_classes.wl`, `verify_extra.wl` | atlas (30 coronas), closure, V(P, n) for the 59 programs of `results/proof_many.json` and 15 further programs, exact class counts |
| `py_extra.py` | the 15 further programs re-checked with the repository's Python code (numpy path) |
| `bigwin.wl`, `survey_exact.wl` | pentagrid windows 46 / 160 / 300 and window surveys with exact congruence classes |
| `figs_basic.wl`, `figs_subst.wl`, `figs_orbits.wl`, `figs_pass.wl`, `figs_tables.wl` | figures and tables of the paper (PNG, 1800 px wide, written to `figs/`) |
| `data/` | small result files (`*.wxf`) and logs of the runs used in the paper |

Run a script with

```
wolfram -noinit -noprompt -script run.wl figs_basic.wl
```

(`run.wl` loads `lib.wl` and the target file as UTF-8; `w.sh` is the same for Git Bash).
`bigwin.wl 300` must be run before `figs_orbits.wl` / `survey_exact.wl` (the cached window is ~170 MB and
is not included).

Notes:

* The numbers of orbit classes written by the GPU survey (`results/surveyall_*.json` etc.) were obtained from a
  floating-point signature (float32) and are too large for programs with big orbits
  (e.g. (SSRSSSSR)*: 38 instead of 13, (SSSR)* on window 300: 63 instead of 41).  The paper uses
  exact congruence classes computed here; closed / escaped counts agree exactly with the GPU survey.
* 15 further programs are proved bounded (single exits k = 35, 40, 45, 50 and 11 exit sequences that
  narrowly missed the candidate threshold of the survey), see `verify_extra.wl`, `py_extra.py`.
