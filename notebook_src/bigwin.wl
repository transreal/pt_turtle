(* bigwin.wl N : build and cache a pentagrid window (same window as survey.py gen N) *)
nwin = ToExpression[$args[[2]]];
t0 = AbsoluteTime[];
rh = ptPentagrid[{-3, -3}, {nwin, nwin}];
Print["window ", nwin, ": ", Length[rh["quad"]], " tiles  [", Round[AbsoluteTime[] - t0], "s]"];
nav = ptNav[rh];
Print["nav built; interior tiles: ", Count[Min /@ Partition[nav["sNext"], 4], _?Positive], "  [", Round[AbsoluteTime[] - t0], "s]"];
xy = ptRhombXY[rh];
Export[FileNameJoin[{$dataDir, "win" <> ToString[nwin] <> ".wxf"}], <|"quad" -> rh["quad"], "fat" -> rh["fat"], "sNext" -> nav["sNext"], "xy" -> Developer`ToPackedArray[xy], "cent" -> Developer`ToPackedArray[Mean /@ xy]|>, PerformanceGoal -> "Speed"];
Print["saved  [", Round[AbsoluteTime[] - t0], "s]"];
