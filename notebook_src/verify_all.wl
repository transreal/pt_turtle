(* independent WL verification of all bounded programs listed in results/proof_many.json *)
keys = Import[FileNameJoin[{$dataDir, "atlas.wxf"}]];
pm = Import[FileNameJoin[{$repoDir, "results", "proof_many.json"}], "RawJSON"];
byN = GroupBy[Keys[pm], pm[#]["n"] &];
Print["programs: ", Length[pm], "  levels: ", KeySort[Length /@ byN]];
out = <||>;
Do[
  Module[{patches, t1 = AbsoluteTime[]},
   patches = Table[Module[{t = ptCanonToTri[k], c = ptCanonCentre[k], st, rh, nav, mask, inside},
      st = ptSubstitute[t, n]; rh = ptRhombs[st]; nav = ptNav[rh];
      mask = Unitize[st["anc"] - c];   (* 0 for descendants of the centre *)
      inside = Flatten[Position[mask[[rh["tri"][[All, 1]]]]*mask[[rh["tri"][[All, 2]]]], 0]];
      {nav["sNext"], Flatten[Table[4 (r - 1) + e + 1, {r, inside}, {e, 0, 3}]]}], {k, keys}];
   Print["level ", n, ": patches built, ", Total[Length[#[[2]]] & /@ patches], " starts  [", Round[AbsoluteTime[] - t1], "s]"];
   Do[Module[{p = ptProg[prog], tot = 0, esc = 0, to = 0, spec = <||>, res},
      Do[res = ptScanC[pt[[1]], p, pt[[2]], 100000];
       tot += Length[pt[[2]]]; esc += Count[res[[All, 1]], 1]; to += Count[res[[All, 1]], 2];
       spec = Merge[{spec, Counts[Select[res, #[[1]] == 0 &][[All, 2]]]}, Total], {pt, patches}];
      spec = KeySort[spec];
      out[prog] = <|"n" -> n, "starts" -> tot, "escape" -> esc, "timeout" -> to, "spectrum" -> spec,
        "agree" -> (tot == pm[prog]["starts"] && esc == 0 && to == 0 && Normal[KeyMap[ToString, spec]] === Normal[KeySortBy[pm[prog]["spectrum"], ToExpression]])|>;
      Print["  (", prog, ")* n=", n, " starts ", tot, " esc ", esc, " to ", to, " agree ", out[prog]["agree"], " periods ", Keys[spec]]], {prog, byN[n]}]],
  {n, Sort[Keys[byN]]}];
Print["ALL AGREE: ", And @@ (#["agree"] & /@ Values[out]), "  count ", Length[out]];
Export[FileNameJoin[{$dataDir, "verify_all.wxf"}], out];
