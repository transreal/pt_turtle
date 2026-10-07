(* verify_extra.wl : V(P, n) for borderline candidates of the exit-sequence survey (third implementation only) *)
keys = Import[FileNameJoin[{$dataDir, "atlas.wxf"}]];
n = ToExpression[$args[[2]]];
progs = Drop[$args, 2];
outFile = FileNameJoin[{$dataDir, "verify_extra.wxf"}];
out = If[FileExistsQ[outFile], Import[outFile], <||>];
patches = Table[Module[{t = ptCanonToTri[k], c = ptCanonCentre[k], st, rh, nav, mask, inside},
    st = ptSubstitute[t, n]; rh = ptRhombs[st]; nav = ptNav[rh];
    mask = Unitize[st["anc"] - c];
    inside = Flatten[Position[mask[[rh["tri"][[All, 1]]]]*mask[[rh["tri"][[All, 2]]]], 0]];
    {nav["sNext"], Flatten[Table[4 (r - 1) + e + 1, {r, inside}, {e, 0, 3}]], rh["quad"]}], {k, keys}];
Print["level ", n, ": patches built"];
Do[Module[{p = ptProg[prog], tot = 0, esc = 0, to = 0, spec = <||>, tr = <||>, cls = <||>, t2 = AbsoluteTime[]},
   Do[Module[{res = ptScanC[pt[[1]], p, pt[[2]], 100000], reps},
     tot += Length[pt[[2]]]; esc += Count[res[[All, 1]], 1]; to += Count[res[[All, 1]], 2];
     spec = Merge[{spec, Counts[Select[res, #[[1]] == 0 &][[All, 2]]]}, Total];
     If[esc == 0 && to == 0,
      reps = DeleteDuplicatesBy[Select[res, #[[1]] == 0 &], #[[3]] &];
      Do[Module[{cells = DeleteDuplicates[ptOrbitRhombsC[pt[[1]], p, r[[3]], r[[2]]]], q, key},
        q = pt[[3]][[cells]]; key = {r[[2]], Hash[ptTransCanon[q]]};
        If[! KeyExistsQ[tr, key], tr[key] = q]], {r, reps}]]], {pt, patches}];
   If[esc == 0 && to == 0, KeyValueMap[Function[{key, q}, cls[{key[[1]], Hash[ptD10Canon[q]]}] = Length[q]], tr]];
   out[prog <> "@" <> ToString[n]] = <|"prog" -> prog, "n" -> n, "starts" -> tot, "escape" -> esc, "timeout" -> to, "spectrum" -> KeySort[spec], "nclasses" -> Length[cls]|>;
   Print["  (", prog, ")* n=", n, " starts ", tot, " esc ", esc, " to ", to, " classes ", Length[cls], " periods ", Keys[KeySort[spec]], "  [", Round[AbsoluteTime[] - t2], "s]"];
   Export[outFile, out]], {prog, progs}];
