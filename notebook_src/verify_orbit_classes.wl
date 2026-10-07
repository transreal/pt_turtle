(* verify_orbit_classes.wl : exact classification of the closed orbits of (P)* found in the verification patches
   sigma^n(C), C in the atlas, as cyclic sequences of states (not just swept sets):
   an orbit is encoded by the sequence of its states at the cycle boundaries, each state by the four vertices of
   its rhomb in Z[zeta] starting at the heading edge (counter-clockwise); two orbits are identified if one is
   mapped onto the other by a translation and a rotation by a multiple of 36 degrees (orientation-preserving
   elements of G) together with a cyclic shift of the sequence.
   usage: run.wl verify_orbit_classes.wl PROG n     (default SSR 7) *)
prog = If[Length[$args] >= 2, $args[[2]], "SSR"]; n = If[Length[$args] >= 3, ToExpression[$args[[3]]], 7];
keys = Import[FileNameJoin[{$dataDir, "atlas.wxf"}]];
p = ptProg[prog]; np = Length[p];
stateQuad[quad_, s_] := With[{q = quad[[Quotient[s - 1, 4] + 1]], e = Mod[s - 1, 4]}, RotateLeft[q, e]];   (* 4 vertices from the heading edge, ccw *)
orbitCode[quad_, sNext_, rep_, per_] := Module[{run = ptRun[sNext, p, rep, per + 2], st},
   st = run["states"]; If[Head[st] === Internal`Bag, st = Internal`BagPart[st, All]];
   st = st[[1 ;; per np ;; np]];           (* states at the cycle boundaries, one period *)
   Developer`ToPackedArray[Table[stateQuad[quad, s], {s, st}]]];                (* per x 4 x 4 integers *)
(* canonical form: for every cyclic shift k, rotate so that the heading edge of the first state points in direction 0
   (a multiple of 36 deg, hence an element of G) and translate its first vertex to the origin; take the lexicographic minimum. *)
unitDirs = Table[ptUnit[j], {j, 0, 9}];
canon[code_] := Module[{per = Length[code], cands},
   cands = Table[Module[{q = code[[k + 1]], r, M, c},
      r = First[FirstPosition[unitDirs, q[[2]] - q[[1]]]] - 1;
      M = ptRhoPow[[Mod[-r, 10] + 1]];
      c = Map[# . M &, RotateLeft[code, k], {2}];
      Flatten[c - ConstantArray[c[[1, 1]], {per, 4}]]], {k, 0, per - 1}];
   First[Sort[cands]]];
tot = 0; esc = 0; to = 0; orbits = <||>; t0 = AbsoluteTime[];
Do[Module[{t = ptCanonToTri[k], c = ptCanonCentre[k], st, rh, nav, mask, inside, starts, res, reps},
   st = ptSubstitute[t, n]; rh = ptRhombs[st]; nav = ptNav[rh];
   mask = Unitize[st["anc"] - c];
   inside = Flatten[Position[mask[[rh["tri"][[All, 1]]]]*mask[[rh["tri"][[All, 2]]]], 0]];
   starts = Flatten[Table[4 (r - 1) + e + 1, {r, inside}, {e, 0, 3}]];
   res = ptScanC[nav["sNext"], p, starts, 100000];
   tot += Length[starts]; esc += Count[res[[All, 1]], 1]; to += Count[res[[All, 1]], 2];
   reps = DeleteDuplicatesBy[Select[res, #[[1]] == 0 &], #[[3]] &];
   Do[Module[{code = orbitCode[rh["quad"], nav["sNext"], r[[3]], r[[2]]], key},
      key = {r[[2]], canon[code]};
      If[KeyExistsQ[orbits, key], orbits[key] = orbits[key] + 1, orbits[key] = 1]], {r, reps}]], {k, keys}];
Print["(", prog, ")* level ", n, ": starts ", tot, " escaped ", esc, " timeout ", to, "  [", Round[AbsoluteTime[] - t0], "s]"];
cls = KeyValueMap[Function[{key, cnt}, <|"period" -> key[[1]], "ncells" -> Length[DeleteDuplicates[key[[2]] ;; ]], "norbits" -> cnt|>], orbits];
Print["distinct closed orbits (up to translation + rotation by multiples of 36 deg, as cyclic state sequences): ", Length[orbits]];
Print["by period: ", KeySort[Counts[Keys[orbits][[All, 1]]]], "   orbits per class: ", SortBy[KeyValueMap[{#1[[1]], #2} &, orbits], First]];
Export[FileNameJoin[{$dataDir, "orbit_classes_" <> prog <> "_n" <> ToString[n] <> ".wxf"}], <|"prog" -> prog, "n" -> n, "starts" -> tot, "escape" -> esc, "timeout" -> to, "classes" -> SortBy[KeyValueMap[{#1[[1]], #2} &, orbits], First]|>];
