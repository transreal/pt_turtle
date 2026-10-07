(* independent WL re-computation: atlas of marked triangle coronas, closure, and verification of bounded programs *)
Print["phi check: ", ptEmbed[{1, 0, 0, 0} . ptPhiM], "  rho: ", ptEmbed[ptUnit[1]]];

atlasFrom[seed_, n_] := Module[{t = ptSubstitute[seed, n], a}, a = ptAtlas[t]; Print["  sigma^", n, ": ", ptNTri[t], " triangles, ", Length[a], " corona classes"]; a];
Print["atlas from wheel:"]; aW = atlasFrom[ptWheel[], 8];
Print["atlas from acute:"]; aA = atlasFrom[ptTri[{ptStdTriangle[0, 1]}], 10];
Print["atlas from obtuse:"]; aO = atlasFrom[ptTri[{ptStdTriangle[1, 1]}], 10];
Print["same key sets: ", Sort[Keys[aW]] === Sort[Keys[aA]], " ", Sort[Keys[aW]] === Sort[Keys[aO]]];
keys = Sort[Keys[aA]];
Print["sizes: ", Sort[Length /@ keys]];
Print["centre types (col, sign): ", Tally[Table[With[{t = ptCanonToTri[k], c = ptCanonCentre[k]}, {t["col"][[c]], ptSign[t][[c]]}], {k, keys}]]];

(* mirror-partner property of the centre triangle in every corona *)
mirrorOK = Table[Module[{t = ptCanonToTri[k], c = ptCanonCentre[k], A, B, C, col, cand},
    {col, A, B, C} = {t["col"][[c]], t["A"][[c]], t["B"][[c]], t["C"][[c]]};
    cand = Select[Range[ptNTri[t]], # != c && Sort[{t["B"][[#]], t["C"][[#]]}] === Sort[{B, C}] &];
    Length[cand] == 1 && t["col"][[First[cand]]] == col && t["A"][[First[cand]]] === B + C - A && t["B"][[First[cand]]] === B && t["C"][[First[cand]]] === C], {k, keys}];
Print["mirror partner property holds in all coronas: ", And @@ mirrorOK];

(* closure *)
keySet = Association[# -> True & /@ keys];
newKeys = {};
Do[Module[{t = ptCanonToTri[k], c = ptCanonCentre[k], st, inc, cor, kk},
   st = ptSubstitute[t, 1]; inc = ptTriIncidence[st];
   Do[If[st["anc"][[i]] == c,
     cor = ptCorona[inc, i];
     If[cor === None, Print["INCOMPLETE corona inside sigma(centre)"]; AppendTo[newKeys, "incomplete"],
      kk = ptCoronaCanon[st, cor, i]; If[! KeyExistsQ[keySet, kk], AppendTo[newKeys, kk]]]], {i, ptNTri[st]}]], {k, keys}];
Print["closure: coronas outside atlas = ", Length[newKeys]];

(* margin: distance from centre triangle to the boundary of its corona *)
Export[FileNameJoin[{$dataDir, "atlas.wxf"}], keys];

(* verification *)
verify[prog_String, n_Integer] := Module[{p = ptProg[prog], tot = 0, esc = 0, to = 0, spec = <||>, per = {}},
   Do[Module[{t = ptCanonToTri[k], c = ptCanonCentre[k], st, rh, nav, inside, starts, res},
     st = ptSubstitute[t, n]; rh = ptRhombs[st]; nav = ptNav[rh];
     inside = Flatten[Position[rh["tri"], {a_, b_} /; (st["anc"][[a]] == c || st["anc"][[b]] == c), {1}, Heads -> False]];
     starts = Flatten[Table[4 (r - 1) + e + 1, {r, inside}, {e, 0, 3}]];
     res = ptScanC[nav["sNext"], p, starts, 200000];
     tot += Length[starts]; esc += Count[res[[All, 1]], 1]; to += Count[res[[All, 1]], 2];
     spec = Merge[{spec, Counts[Select[res, #[[1]] == 0 &][[All, 2]]]}, Total];
     AppendTo[per, {Length[k], Length[rh["quad"]], Length[starts], Count[res[[All, 1]], 0]}]], {k, keys}];
   Print["(", prog, ")* n=", n, ": starts ", tot, " escapes ", esc, " timeouts ", to, " spectrum ", KeySort[spec]];
   <|"prog" -> prog, "n" -> n, "starts" -> tot, "escape" -> esc, "timeout" -> to, "spectrum" -> KeySort[spec], "per" -> per|>];
r1 = verify["SSR", 7];
r2 = verify["SR", 7];
r3 = verify["SRSRSRSSR", 7];
r4 = verify["SSRSRSRSRSSR", 7];
Export[FileNameJoin[{$dataDir, "verify_n7.wxf"}], {r1, r2, r3, r4}];
