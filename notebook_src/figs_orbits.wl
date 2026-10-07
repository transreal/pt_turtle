(* figs_orbits.wl : orbit-class galleries, hierarchical rings, diffusive trajectory (window 300) *)
tx[s_, pos_, size_ : 9, opts___] := Text[ptText[s, size, opts], pos];
title[s_] := ptText[s, 9.5, Black];
win = Import[FileNameJoin[{$dataDir, "win300.wxf"}]];
sNext = win["sNext"]; cent = win["cent"]; quad = win["quad"]; xy = win["xy"]; rhF = <|"fat" -> win["fat"]|>;
interior = Flatten[Position[Min /@ Partition[sNext, 4], _?Positive]];
starts = Flatten[Outer[Plus, 4 (interior - 1), {1, 2, 3, 4}]];
ctr = Mean[cent];
exact[prog_String] := Module[{f = FileNameJoin[{$dataDir, "exact300_" <> prog <> ".wxf"}], res, cls},
   If[FileExistsQ[f], Import[f],
    res = ptScanC[sNext, ptProg[prog], starts, 100000]; cls = ptExactClasses[quad, sNext, cent, ptProg[prog], res]; Export[f, cls]; cls]];
cellsOf[prog_String, c_] := DeleteDuplicates[ptOrbitRhombsC[sNext, ptProg[prog], c["rep"], c["period"]]];
gallery[name_, prog_String, ncol_, size_ : 600] := Module[{cls = exact[prog], panels},
   Print[name, ": (", prog, ")* classes ", Length[cls], "  ", {#period, #ncells} & /@ cls];
   panels = Table[ptPanel[rhF, xy, cent, cellsOf[prog, c], "Margin" -> 0.8, "Label" -> "(" <> ToString[c["period"]] <> ", " <> ToString[c["ncells"]] <> ")", "LabelSize" -> 8.5, "Edges" -> (c["ncells"] < 200)], {c, cls}];
   panels = Partition[panels, UpTo[ncol]]; panels[[-1]] = Join[panels[[-1]], ConstantArray[Graphics[{}], ncol - Length[panels[[-1]]]]];
   ptSave[name, GraphicsGrid[panels, Spacings -> {4, 6}, ImageSize -> size]]];
gallery["classes-ssr", "SSR", 4];
gallery["classes-6", "SRSRSRSSR", 4];
gallery["classes-36", "SSRSRSRSRSSR", 4];
gallery["classes-3322", "SSRSSSSR", 4];

(* ---------- (SSSR)*: concentric families ---------- *)
clsS = exact["SSSR"];
resS = ptScanC[sNext, ptProg["SSSR"], starts, 100000];
orbS = ptOrbitData[sNext, cent, ptProg["SSSR"], resS, 40];
Print["SSSR orbits with period >= 40: ", Length[orbS]];
big = First[SortBy[Select[orbS, #period == 1385 &], Norm[#centre - ctr] &]];
bundle = SortBy[Select[orbS, Norm[#centre - big["centre"]] < 0.01 &], #period &];
Print["bundle around the period-1385 orbit nearest the window centre: periods ", #period & /@ bundle, " radii ", Round[#radius & /@ bundle, 0.1]];
Print["all bundles (by centre) with >= 4 orbits: ", Take[Reverse[SortBy[Tally[Round[#centre, 0.01] & /@ orbS], Last]], UpTo[5]]];
Export[FileNameJoin[{$dataDir, "sssr_bundle.wxf"}], bundle];
ringFig[sel_, size_, lab_ : True] := Module[{cols = {RGBColor[0.85, 0.15, 0.15], RGBColor[0.95, 0.55, 0.1], RGBColor[0.15, 0.6, 0.3], RGBColor[0.1, 0.45, 0.85], RGBColor[0.55, 0.25, 0.75], RGBColor[0.2, 0.2, 0.2], RGBColor[0.7, 0.5, 0.2], RGBColor[0.1, 0.65, 0.7]}, cellsL, R},
   cellsL = cellsOf["SSSR", #] & /@ sel; R = Max[#radius & /@ sel] + 2;
   Graphics[{EdgeForm[None], MapIndexed[{cols[[Mod[#2[[1]] - 1, 8] + 1]], Polygon[xy[[#1]]]} &, cellsL]},
    PlotRange -> {big["centre"][[1]] + {-R, R}, big["centre"][[2]] + {-R, R}}, ImageSize -> size]];
ringCols3 = {RGBColor[0.15, 0.6, 0.3], RGBColor[0.95, 0.55, 0.1], RGBColor[0.85, 0.15, 0.15]};
ringsSel = Select[bundle, MemberQ[{85, 345, 1385, 5545}, #period] &];
Print["rings selected: ", {#period, #ncells, Round[#radius, 0.1]} & /@ ringsSel];
pathFig[sel_, size_] := Module[{cols = {RGBColor[0.85, 0.15, 0.15], RGBColor[0.1, 0.45, 0.85], RGBColor[0.15, 0.6, 0.3], RGBColor[0.95, 0.55, 0.1]}, R},
   R = Max[#radius & /@ sel] + 1.5;
   Graphics[MapIndexed[{cols[[#2[[1]]]], AbsoluteThickness[0.6], Line[Append[#, First[#]] &[cent[[ptOrbitRhombsC[sNext, ptProg["SSSR"], #1["rep"], #1["period"]]]]]]} &, sel],
    PlotRange -> {big["centre"][[1]] + {-R, R}, big["centre"][[2]] + {-R, R}}, ImageSize -> size]];
sameScale = Select[bundle, MemberQ[{290, 345, 380, 615}, #period] &];
Print["same-scale orbits: ", {#period, #ncells, Round[#radius, 0.1]} & /@ sameScale];
ptSave["sssr-rings", GraphicsRow[{Show[ringFig[Reverse[ringsSel], Automatic], PlotLabel -> title["(a) 周期 85, 345, 1385, 5545 の軌道"]], Show[pathFig[sameScale, Automatic], PlotLabel -> title["(b) 周期 290, 345, 380, 615 の軌道 (同じ中心)"]]}, Spacings -> 8, ImageSize -> 600]];

(* families of periods *)
Module[{fam = {{45, 85, 175, 345, 695, 1385, 2775, 5545}, {40, 100, 180, 380, 740, 1500, 2980, 5980}, {20, 110, 290, 510, 1090, 2110, 4290}, {265, 85, 615, 785, 2015, 3585, 7615}}, have = #period & /@ clsS},
  Print["families present in window 300: ", Map[MemberQ[have, #] &, fam, {2}]];
  ptSave["sssr-families", ListLogPlot[MapIndexed[Function[{f, i}, Transpose[{Range[Length[f]], f}]], fam], Joined -> True, PlotMarkers -> {Automatic, 5}, Frame -> True,
    FrameLabel -> {ptText["列の中での順番", 9], ptText["周期 (サイクル数)", 9]}, PlotLegends -> Placed[LineLegend[{"2P ± 5", "2P ± 20", "2P ± 70", "2P ± 445"}, LabelStyle -> {FontFamily -> ptFont, FontSize -> 8}], {0.22, 0.78}],
    Epilog -> {GrayLevel[0.3], Dashed, Line[{{1, Log[12.]}, {8, Log[12. 2^7]}}]}, ImageSize -> 430, LabelStyle -> {FontFamily -> ptFont, FontSize -> 8, Black}, GridLines -> Automatic, GridLinesStyle -> GrayLevel[0.9], PlotRange -> All]]];

(* ---------- diffusive trajectory (SSLSR)* from the centre ---------- *)
c0 = First[Ordering[Norm /@ ((# - ctr) & /@ cent[[interior]]), 1]]; c0 = interior[[c0]];
trajs = Table[ptRun[sNext, ptProg["SSLSR"], 4 (c0 - 1) + e + 1, 200000], {e, 0, 3}];
Print["(SSLSR)* from the centre tile ", c0, " at ", cent[[c0]], ": {flag, cycles, tiles} = ", {#flag, #cycles, Length[Union[#rhombs]]} & /@ trajs];
tr = First[SortBy[Select[trajs, #flag == 1 &], -#cycles &]];
Print["chosen heading edge: ", First[FirstPosition[trajs, tr]] - 1];
trajFig[size_] := Module[{cells = Union[tr["rhombs"]], cyc},
   cyc = cent[[ptStateRhomb /@ tr["states"][[1 ;; -1 ;; 5]]]];
   Graphics[{{EdgeForm[None], ptBlue, Polygon[xy[[cells]]]}, {ptRed, AbsoluteThickness[0.25], Line[cyc]}, {Black, PointSize[0.012], Point[cent[[c0]]]},
     {GrayLevel[0.4], AbsoluteThickness[0.6], Line[{{-7.5, -7.5}, {742.5, -7.5}, {742.5, 742.5}, {-7.5, 742.5}, {-7.5, -7.5}}]}}, PlotRange -> {{-20, 755}, {-20, 755}}, ImageSize -> size]];
Print["trajectory used: cycles ", tr["cycles"], " tiles ", Length[Union[tr["rhombs"]]], " window bbox ", MinMax[cent[[All, 1]]], MinMax[cent[[All, 2]]]];
ptSave["diffusive-trajectory", Show[trajFig[480], ImageSize -> 480]];
Export[FileNameJoin[{$dataDir, "traj_info.wxf"}], <|"cycles" -> tr["cycles"], "tiles" -> Length[Union[tr["rhombs"]]]|>];

(* ---------- rms displacement ---------- *)
depth = With[{x = cent[[All, 1]], y = cent[[All, 2]]}, MapThread[Min, {x - Min[x], Max[x] - x, y - Min[y], Max[y] - y}]];
central = Flatten[Position[depth, _?(# >= 145 &)]]; central = Intersection[central, interior];
SeedRandom[20261005];
chk = {50, 100, 200, 400, 800, 1600, 3200, 6400};
msd[prog_String, ntr_ : 200] := Module[{p = ptProg[prog], st = RandomChoice[central, ntr], ed = RandomInteger[{0, 3}, ntr], runs, D, ncl = 0, nesc = 0, rms, fit},
   runs = MapThread[Function[{c, e}, {c, ptRunCheckC[sNext, p, 4 (c - 1) + e + 1, chk]}], {st, ed}];
   ncl = Count[runs, {_, r_} /; r[[-2]] == 0]; nesc = Count[runs, {_, r_} /; r[[-2]] == 1];
   D = Table[Select[Table[If[r[[2, -2]] != 0 && r[[2, k]] > 0, Norm[cent[[ptStateRhomb[r[[2, k]]]]] - cent[[r[[1]]]]], Nothing], {r, runs}], NumericQ], {k, Length[chk]}];
   rms = Table[If[Length[D[[k]]] >= 5, {chk[[k]], Sqrt[Mean[D[[k]]^2]]}, Nothing], {k, Length[chk]}];
   fit = Fit[Log[rms], {1, x}, x];
   Print["(", prog, ")*: closed ", ncl, " escaped ", nesc, " of ", ntr, "; rms ", Round[rms, 0.1], "  n ", Length /@ D, "  exponent ", Coefficient[fit, x]];
   <|"prog" -> prog, "rms" -> rms, "exponent" -> Coefficient[fit, x], "closed" -> ncl, "escaped" -> nesc|>];
ms = msd /@ {"SSLSR", "SSSSLSR", "SSRSSSSSSR"};
Export[FileNameJoin[{$dataDir, "msd.wxf"}], ms];
ptSave["diffusion-msd", ListLogLogPlot[#rms & /@ ms, Joined -> True, PlotMarkers -> {Automatic, 6}, Frame -> True, FrameLabel -> {ptText["サイクル数 t", 9], ptText["変位の二乗平均平方根", 9]},
   PlotLegends -> Placed[LineLegend[{"(SSLSR)*", "(SSSSLSR)*", "(SSR)(S⁶R)"}, LabelStyle -> {FontFamily -> ptFont, FontSize -> 8}], {0.25, 0.8}],
   Epilog -> {GrayLevel[0.3], Dashed, Line[{{Log[50.], Log[12.]}, {Log[6400.], Log[12. Sqrt[128.]]}}], Text[ptText["傾き 1/2", 8], {Log[1500.], Log[40.]}]},
   ImageSize -> 430, LabelStyle -> {FontFamily -> ptFont, FontSize -> 8, Black}, GridLines -> Automatic, GridLinesStyle -> GrayLevel[0.9], PlotRange -> All]];

(* ---------- teaser: three regimes ---------- *)
Module[{cS = Last[exact["SSR"]], g1, g2, g3},
  g1 = ptPanel[rhF, xy, cent, cellsOf["SSR", cS], "Margin" -> 3.5, "Label" -> "(a) 一様有界: (SSR)*", "LabelSize" -> 9];
  g2 = Show[ringFig[Reverse[Select[bundle, MemberQ[{85, 345, 1385}, #period] &]], Automatic], PlotLabel -> ptText["(b) 階層型: (SSSR)*", 9, Black]];
  g3 = Show[trajFig[Automatic], PlotLabel -> ptText["(c) 拡散型: (SSLSR)*", 9, Black]];
  ptSave["three-regimes", GraphicsRow[{g1, g2, g3}, Spacings -> 6, ImageSize -> 600]]];
