(* numbers requested by the review: initial states per window, closed fractions of (SSSR)*, corona margin delta0 *)
Do[Module[{win = Import[FileNameJoin[{$dataDir, "win" <> ToString[w] <> ".wxf"}]], interior},
   interior = Flatten[Position[Min /@ Partition[win["sNext"], 4], _?Positive]];
   Print["window ", w, ": tiles ", Length[win["cent"]], "  start tiles (4 neighbours inside) ", Length[interior], "  initial states ", 4 Length[interior]]], {w, {46, 160, 300}}];
Do[Module[{d = Import[FileNameJoin[{$dataDir, "survey_exact_all43_" <> ToString[w] <> ".wxf"}]], r},
   r = SelectFirst[d, (#["prog"] === "SSSR" || Lookup[#, "program", ""] === "SSSR") &];
   If[AssociationQ[r], Print["window ", w, " (SSSR)*: ", KeyTake[r, {"prog", "program", "closed", "escaped", "timeout", "starts", "nstates", "closedFrac", "nclasses", "maxPeriod", "maxP"}]],
     Print["window ", w, ": keys of first record ", Keys[First[d]]]]], {w, {46, 160, 300}}];
keys = Import[FileNameJoin[{$dataDir, "atlas.wxf"}]];
segDist[p_, {a_, b_}] := Module[{d = b - a, t}, t = Clip[(p - a) . d/(d . d), {0, 1}]; Norm[p - (a + t d)]];
delta0 = Table[Module[{t = ptCanonToTri[k], c = ptCanonCentre[k], xy, edges, bnd, cv},
     xy = ptTriXY[t]; edges = Flatten[Table[Sort /@ Partition[Append[Round[xy[[j]], 10.^-8], Round[xy[[j, 1]], 10.^-8]], 2, 1], {j, Length[xy]}], 1];
     bnd = Select[Tally[edges], #[[2]] == 1 &][[All, 1]]; cv = xy[[c]];
     Min[Join[Flatten[Table[segDist[p, e], {p, cv}, {e, bnd}]], Flatten[Table[segDist[q, {cv[[i]], cv[[Mod[i, 3] + 1]]}], {q, Flatten[bnd, 1]}, {i, 3}]]]]], {k, keys}];
Print["corona margins: min ", Min[delta0], "  max ", Max[delta0], "  sin 36 = ", Sin[36. Degree], "  all >= sin36 - 1e-9: ", And @@ Thread[delta0 >= Sin[36. Degree] - 10.^-9]];
