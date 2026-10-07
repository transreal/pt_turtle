(* survey_exact.wl : independent re-computation of the window surveys with exact orbit classes.
   usage: survey_exact.wl N set      set = all43 | single | pairs *)
nwin = ToExpression[$args[[2]]]; which = $args[[3]];
f = FileNameJoin[{$dataDir, "win" <> ToString[nwin] <> ".wxf"}];
If[! FileExistsQ[f],
  Module[{rh = ptPentagrid[{-3, -3}, {nwin, nwin}], nav, xy}, nav = ptNav[rh]; xy = ptRhombXY[rh];
   Export[f, <|"quad" -> rh["quad"], "fat" -> rh["fat"], "sNext" -> nav["sNext"], "xy" -> Developer`ToPackedArray[xy], "cent" -> Developer`ToPackedArray[Mean /@ xy]|>, PerformanceGoal -> "Speed"]]];
win = Import[f]; sNext = win["sNext"]; cent = win["cent"]; quad = win["quad"];
interior = Flatten[Position[Min /@ Partition[sNext, 4], _?Positive]];
starts = Flatten[Outer[Plus, 4 (interior - 1), {1, 2, 3, 4}]];
Print["window ", nwin, ": tiles ", Length[cent], " starts ", Length[starts]];
single[k_] := StringJoin[ConstantArray["SR", k - 3]] <> "SSR";
progs = Switch[which,
   "all43", #["prog"] & /@ Import[FileNameJoin[{$repoDir, "results", "surveyall_300_8.json"}], "RawJSON"],
   "single", single /@ Join[Range[3, 30], {35, 40, 45, 50}],
   "pairs", #["prog"] & /@ Select[Import[FileNameJoin[{$repoDir, "results", "surveyexits_300.json"}], "RawJSON"], Length[#["exits"]] == 2 &]];
ref = If[which === "all43", Association[#["prog"] -> # & /@ Import[FileNameJoin[{$repoDir, "results", "surveyall_" <> ToString[nwin] <> "_8.json"}], "RawJSON"]], <||>];
outFile = FileNameJoin[{$dataDir, "survey_exact_" <> which <> "_" <> ToString[nwin] <> ".wxf"}];
out = If[FileExistsQ[outFile], Import[outFile], <||>];
Do[If[! KeyExistsQ[out, prog],
   Module[{p = ptProg[prog], res, cls, t0 = AbsoluteTime[], rec},
    res = ptScanC[sNext, p, starts, 100000];
    cls = ptExactClasses[quad, sNext, cent, p, res];
    rec = <|"prog" -> prog, "N" -> nwin, "starts" -> Length[starts], "closed" -> Count[res[[All, 1]], 0], "escape" -> Count[res[[All, 1]], 1], "timeout" -> Count[res[[All, 1]], 2],
      "nclasses" -> Length[cls], "ntrans" -> Total[#ntrans & /@ cls], "norbits" -> Total[#norbits & /@ cls], "maxperiod" -> If[cls === {}, 0, Max[#period & /@ cls]], "maxradius" -> If[cls === {}, 0., Max[#radius & /@ cls]],
      "classes" -> ({#period, #ncells, #norbits} & /@ cls)|>;
    out[prog] = rec;
    Print["(", prog, ")* closed ", rec["closed"], " esc ", rec["escape"], " to ", rec["timeout"], " classes ", rec["nclasses"], " maxP ", rec["maxperiod"], " maxR ", Round[rec["maxradius"], 0.1],
     If[KeyExistsQ[ref, prog], "  | ref closed " <> ToString[ref[prog]["closed"]] <> If[ref[prog]["closed"] == rec["closed"], " OK", " MISMATCH"] <> " ref classes " <> ToString[ref[prog]["nclasses"]], ""], "  [", Round[AbsoluteTime[] - t0], "s]"];
    Export[outFile, out]]], {prog, progs}];
Print["finished ", which, " ", nwin];
