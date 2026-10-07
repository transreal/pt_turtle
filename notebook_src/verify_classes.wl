(* verify_classes.wl : exact congruence classes of the orbits found in the verification patches (all bounded programs),
   plus the extra single-exit programs k = 35, 40, 45, 50 and checks (b) of the atlas *)
keys = Import[FileNameJoin[{$dataDir, "atlas.wxf"}]];
(* condition (b): every edge of the centre is shared, as a whole edge, with exactly one other triangle of the corona *)
edgeOK = Table[Module[{t = ptCanonToTri[k], c = ptCanonCentre[k], crn, others},
    crn = {t["A"][[c]], t["B"][[c]], t["C"][[c]]};
    others = Table[{t["A"][[j]], t["B"][[j]], t["C"][[j]]}, {j, Complement[Range[ptNTri[t]], {c}]}];
    And @@ Table[Count[others, o_ /; SubsetQ[o, e]] == 1, {e, Subsets[crn, {2}]}]], {k, keys}];
Print["condition (b) holds for all 30 centres: ", And @@ edgeOK];
(* margin delta0 : distance from the centre triangle to the boundary of the corona *)
segDist[p_, {a_, b_}] := Module[{d = b - a, t}, t = Clip[(p - a) . d/(d . d), {0, 1}]; Norm[p - (a + t d)]];
delta0 = Min[Table[Module[{t = ptCanonToTri[k], c = ptCanonCentre[k], xy, edges, bnd, cv},
     xy = ptTriXY[t]; edges = Flatten[Table[Sort /@ Partition[Append[Round[xy[[j]], 10.^-8], Round[xy[[j, 1]], 10.^-8]], 2, 1], {j, Length[xy]}], 1];
     bnd = Select[Tally[edges], #[[2]] == 1 &][[All, 1]]; cv = xy[[c]];
     Min[Join[Flatten[Table[segDist[p, e], {p, cv}, {e, bnd}]], Flatten[Table[segDist[q, {cv[[i]], cv[[Mod[i, 3] + 1]]}], {q, Flatten[bnd, 1]}, {i, 3}]]]]], {k, keys}]];
Print["delta0 = ", delta0, "  sin 36 = ", Sin[36. Degree]];

pm = Import[FileNameJoin[{$repoDir, "results", "proof_many.json"}], "RawJSON"];
extra = Association[Table[StringJoin[ConstantArray["SR", k - 3]] <> "SSR" -> <|"n" -> 7|>, {k, {35, 40, 45, 50}}]];
all = Join[pm, extra];
byN = GroupBy[Keys[all], all[#]["n"] &];
outFile = FileNameJoin[{$dataDir, "verify_classes.wxf"}];
out = If[FileExistsQ[outFile], Import[outFile], <||>];
Do[
  If[! SubsetQ[Keys[out], byN[n]],
   Module[{patches, t1 = AbsoluteTime[]},
    patches = Table[Module[{t = ptCanonToTri[k], c = ptCanonCentre[k], st, rh, nav, mask, inside},
       st = ptSubstitute[t, n]; rh = ptRhombs[st]; nav = ptNav[rh];
       mask = Unitize[st["anc"] - c];
       inside = Flatten[Position[mask[[rh["tri"][[All, 1]]]]*mask[[rh["tri"][[All, 2]]]], 0]];
       {nav["sNext"], Flatten[Table[4 (r - 1) + e + 1, {r, inside}, {e, 0, 3}]], rh["quad"]}], {k, keys}];
    Print["level ", n, ": patches built  [", Round[AbsoluteTime[] - t1], "s]"];
    Do[If[! KeyExistsQ[out, prog],
      Module[{p = ptProg[prog], tot = 0, esc = 0, to = 0, spec = <||>, tr = <||>, cls = <||>, t2 = AbsoluteTime[]},
       Do[Module[{res = ptScanC[pt[[1]], p, pt[[2]], 100000], reps},
         tot += Length[pt[[2]]]; esc += Count[res[[All, 1]], 1]; to += Count[res[[All, 1]], 2];
         spec = Merge[{spec, Counts[Select[res, #[[1]] == 0 &][[All, 2]]]}, Total];
         reps = DeleteDuplicatesBy[Select[res, #[[1]] == 0 &], #[[3]] &];
         Do[Module[{cells = DeleteDuplicates[ptOrbitRhombsC[pt[[1]], p, r[[3]], r[[2]]]], q, key},
           q = pt[[3]][[cells]]; key = {r[[2]], Hash[ptTransCanon[q]]};
           If[! KeyExistsQ[tr, key], tr[key] = q]], {r, reps}]], {pt, patches}];
       KeyValueMap[Function[{key, q}, cls[{key[[1]], Hash[ptD10Canon[q]]}] = Length[q]], tr];
       out[prog] = <|"n" -> n, "starts" -> tot, "escape" -> esc, "timeout" -> to, "spectrum" -> KeySort[spec], "nclasses" -> Length[cls], "ntrans" -> Length[tr],
         "classes" -> Sort[KeyValueMap[{#1[[1]], #2} &, cls]]|>;
       Print["  (", prog, ")* n=", n, " starts ", tot, " esc ", esc, " to ", to, " classes ", Length[cls], " trans ", Length[tr], " periods ", Keys[KeySort[spec]], "  [", Round[AbsoluteTime[] - t2], "s]"];
       Export[outFile, out]]], {prog, byN[n]}]]],
  {n, Sort[Keys[byN]]}];
Print["done: ", Length[out], " programs; all closed: ", And @@ (#["escape"] == 0 && #["timeout"] == 0 & /@ Values[out])];
