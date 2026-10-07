(* figs_basic.wl : definitions, turtle commands, ribbons, pivot walk.
   All figures are composed at ImageSize 600 and rasterised to 1800 px. *)
rh = ptPentagrid[{-9, -9}, {18, 18}];
xy = ptRhombXY[rh]; cent = Mean /@ xy; nav = ptNav[rh]; nR = Length[xy];
near[p_, r_] := Flatten[Position[Norm /@ ((# - p) & /@ cent), _?(# < r &)]];
nearest[p_] := First[Ordering[Norm /@ ((# - p) & /@ cent), 1]];
bgTiles[idx_] := {EdgeForm[{GrayLevel[0.6], AbsoluteThickness[0.4]}], {ptColFat, Polygon[xy[[Select[idx, rh["fat"][[#]] == 1 &]]]]}, {ptColThin, Polygon[xy[[Select[idx, rh["fat"][[#]] == 0 &]]]]}};
st[r_, e_] := 4 (r - 1) + e + 1;
edgeXY[s_] := ptStateDExy[rh, s];
mid[s_] := Mean[edgeXY[s]];
tx[s_, pos_, size_ : 9, opts___] := Text[ptText[s, size, opts], pos];
mx[s_, pos_, size_ : 11, opts___] := Text[ptMath[s, size, opts], pos];
title[s_] := ptText[s, 9.5, Black];
hlTile = RGBColor[0.75, 0.87, 1.];
(* glyph of a turtle state: heading edge + arrow from the centroid *)
glyph[s_, col_] := With[{c = cent[[ptStateRhomb[s]]], m = mid[s]}, {col, AbsoluteThickness[2.6], CapForm["Round"], Line[edgeXY[s]], AbsoluteThickness[1.3], Arrowheads[0.06], Arrow[{c - 0.15 (m - c), c + 0.72 (m - c)}]}];
frame[g_, c_, r_, opts___] := Graphics[g, PlotRange -> {c[[1]] + {-r, r}, c[[2]] + {-r, r}}, PlotRangeClipping -> True, PlotRangePadding -> 0, opts];
frameW[g_, c_, rx_, ry_, opts___] := Graphics[g, PlotRange -> {c[[1]] + {-rx, rx}, c[[2]] + {-ry, ry}}, PlotRangeClipping -> True, PlotRangePadding -> 0, opts];

(* ---------- pentagrid and dual tiling with two ribbons ---------- *)
Module[{w = 2.3, rhs, xys, famCol, lines, hl = {{1, 0}, {3, 1}}, g1, g2, isR},
  rhs = ptPentagrid[{-w, -w}, {2 w, 2 w}]; xys = ptRhombXY[rhs];
  famCol = {RGBColor[0.85, 0.25, 0.2], RGBColor[0.2, 0.5, 0.85], RGBColor[0.2, 0.62, 0.3], RGBColor[0.85, 0.55, 0.1], RGBColor[0.55, 0.3, 0.75]};
  lines = Flatten[Table[With[{d = m + ptOffsets[[j]], e = ptEV5[[j]], p = {-ptEV5[[j, 2]], ptEV5[[j, 1]]}},
      {famCol[[j]], If[MemberQ[hl, {j, m}], AbsoluteThickness[2.2], AbsoluteThickness[0.6]], Line[{d e - 10 p, d e + 10 p}]}], {j, 5}, {m, -5, 5}], 1];
  g1 = Graphics[{lines}, PlotRange -> {{-w, w}, {-w, w}}, PlotRangeClipping -> True, Frame -> True, FrameTicks -> None, FrameStyle -> GrayLevel[0.4]];
  isR[k_] := Function[l, (l[[1]] == hl[[k, 1]] && l[[2]] == hl[[k, 2]]) || (l[[3]] == hl[[k, 1]] && l[[4]] == hl[[k, 2]])] /@ rhs["lines"];
  g2 = Graphics[{EdgeForm[{GrayLevel[0.45], AbsoluteThickness[0.4]}],
     {ptColFat, Polygon[xys[[Flatten[Position[rhs["fat"], 1]]]]]}, {ptColThin, Polygon[xys[[Flatten[Position[rhs["fat"], 0]]]]]},
     {Opacity[0.55], famCol[[hl[[1, 1]]]], Polygon[xys[[Flatten[Position[isR[1], True]]]]]},
     {Opacity[0.55], famCol[[hl[[2, 1]]]], Polygon[xys[[Flatten[Position[isR[2], True]]]]]}}, PlotRange -> {{-7, 7}, {-7, 7}}, PlotRangeClipping -> True, Frame -> True, FrameTicks -> None, FrameStyle -> GrayLevel[0.4]];
  ptSave["pentagrid-ribbons", GraphicsRow[{g1, g2}, Spacings -> 8, ImageSize -> 600]]];

c0 = nearest[{0.6, 0.4}]; e0 = 0;

(* ---------- turtle state and the three commands ---------- *)
Module[{s = st[c0, e0], sS, sR, sL, R = 1.75, ctr, loc, panel},
  sS = ptApply[nav, s, 0]; sR = ptApply[nav, s, 1]; sL = ptApply[nav, s, 2];
  ctr = Mean[{cent[[c0]], cent[[ptStateRhomb[sS]]]}] + {0, 0.15};
  loc = near[ctr, R + 1.8];
  panel[ttl_, hlTiles_, items_] := frame[{bgTiles[loc], {hlTile, EdgeForm[{GrayLevel[0.4], AbsoluteThickness[0.5]}], Polygon[xy[[hlTiles]]]}, items}, ctr, R, PlotLabel -> title[ttl]];
  ptSave["turtle-commands", GraphicsGrid[{{
      panel["(a) 状態 (c, e)", {c0}, {glyph[s, ptRed],
        Table[tx[ToString[k], mid[st[c0, k]] + 0.2 (cent[[c0]] - mid[st[c0, k]]), 8.5, GrayLevel[0.1]], {k, 0, 3}],
        Table[{Black, PointSize[0.022], Point[xy[[c0, k + 1]]], Text[Row[{ptMath["v", 10], Subscript["", ptText[ToString[k], 7]]}], xy[[c0, k + 1]] + 0.2 Normalize[xy[[c0, k + 1]] - cent[[c0]]]]}, {k, 0, 3}]}],
      panel["(b) S: 辺を渡り対辺を新しい進行辺にする", {c0, ptStateRhomb[sS]}, {glyph[s, GrayLevel[0.55]], glyph[sS, ptRed]}]}, {
      panel["(c) R: 進行辺を時計回りに 1 つ回す", {c0}, {glyph[s, GrayLevel[0.55]], glyph[sR, ptRed]}],
      panel["(d) L: 進行辺を反時計回りに 1 つ回す", {c0}, {glyph[s, GrayLevel[0.55]], glyph[sL, ptRed]}]}}, Spacings -> {6, 6}, ImageSize -> 560]]];

(* ---------- ribbon ---------- *)
ribbon[s_, nf_, nb_] := Module[{f = NestWhileList[ptApply[nav, #, 0] &, s, # != 0 &, 1, nf], b, inv},
   inv[x_] := Module[{y = ptApply[nav, ptApply[nav, x, 1], 1]}, y = ptApply[nav, y, 0]; If[y == 0, 0, ptApply[nav, ptApply[nav, y, 1], 1]]];
   b = NestWhileList[inv, s, # != 0 &, 1, nb];
   DeleteCases[Join[Reverse[Rest[b]], f], 0]];
Module[{s = st[c0, 1], rb, tiles, R = 6.2, ctr, loc, d, nrm, s2, rb2},
  rb = ribbon[s, 12, 12]; tiles = ptStateRhomb /@ rb;
  s2 = st[c0, 0]; rb2 = ribbon[s2, 12, 12];
  ctr = cent[[c0]]; loc = near[ctr, R + 3];
  d = Normalize[Subtract @@ Reverse[edgeXY[s]]]; nrm = {d[[2]], -d[[1]]};
  If[nrm . (cent[[ptStateRhomb[ptApply[nav, s, 0]]]] - cent[[c0]]) < 0, nrm = -nrm];
  ptSave["ribbon", frameW[{bgTiles[loc],
      {RGBColor[0.75, 0.9, 0.78], Opacity[0.9], EdgeForm[{GrayLevel[0.45], AbsoluteThickness[0.4]}], Polygon[xy[[ptStateRhomb /@ rb2]]]},
      {RGBColor[1., 0.85, 0.6], EdgeForm[{GrayLevel[0.45], AbsoluteThickness[0.4]}], Polygon[xy[[tiles]]]},
      {RGBColor[0.75, 0.3, 0.1], AbsoluteThickness[2.2], CapForm["Round"], Line[edgeXY /@ rb]},
      {RGBColor[0.1, 0.45, 0.2], AbsoluteThickness[2.2], CapForm["Round"], Line[edgeXY /@ rb2]},
      {ptBlue, AbsoluteThickness[1.1], Arrowheads[{{0.018, 0.6}}], Arrow /@ Partition[cent[[tiles]], 2, 1]},
      {Black, AbsoluteThickness[1.4], Arrowheads[0.03], Arrow[{ctr - 4.6 d - 4.2 nrm, ctr - 4.6 d - 1.8 nrm}], mx["h", ctr - 4.6 d - 1.5 nrm, 13]}}, ctr, R, 3.6, ImageSize -> 600]]];

(* ---------- pivot rules ---------- *)
Module[{s = st[c0, e0], R = 1.65, ctr, loc, panel, dedge},
  ctr = cent[[c0]] + 0.3 (cent[[ptStateRhomb[ptApply[nav, s, 0]]]] - cent[[c0]]); loc = near[ctr, R + 1.8];
  dedge[x_, col_] := {col, AbsoluteThickness[2.4], Arrowheads[0.085], Arrow[edgeXY[x], {0.0, 0.07}], PointSize[0.04], Point[First[edgeXY[x]]]};
  panel[ttl_, c_, extraTiles_] := Module[{t = ptApply[nav, s, c], p = First[edgeXY[s]], p2},
    p2 = First[edgeXY[t]];
    frame[{bgTiles[loc], {hlTile, EdgeForm[{GrayLevel[0.4], AbsoluteThickness[0.5]}], Polygon[xy[[Join[{c0}, extraTiles]]]]},
      dedge[s, RGBColor[0.1, 0.25, 0.7]], dedge[t, ptRed],
      {ptGreen, Dashed, AbsoluteThickness[1.4], Arrowheads[0.06], Arrow[{p + 0.2 Normalize[cent[[c0]] - p], p2 + 0.2 Normalize[cent[[ptStateRhomb[t]]] - p2]}]},
      mx["p", p + 0.26 Normalize[p - ctr - {0, 0.3}], 12, RGBColor[0.1, 0.25, 0.7]], mx["p'", p2 + 0.26 Normalize[p2 - ctr + {0.2, 0.1}], 12, ptRed],
      mx["c", cent[[c0]] + {0.05, 0.22}, 12], If[extraTiles =!= {}, mx["c'", cent[[First[extraTiles]]] + {0.25, 0.}, 12], {}]},
     ctr, R, PlotLabel -> title[ttl]]];
  ptSave["pivot-rules", GraphicsRow[{panel["S", 0, {ptStateRhomb[ptApply[nav, s, 0]]}], panel["L", 2, {}], panel["R", 1, {}]}, Spacings -> 6, ImageSize -> 600]]];

(* ---------- vertex data for stars ---------- *)
vid = nav["vid"];
vpos = Association[Flatten[Table[vid[[r, k]] -> xy[[r, k]], {r, nR}, {k, 4}]]];
vnb = Merge[Flatten[Table[{vid[[r, k]] -> vid[[r, Mod[k, 4] + 1]], vid[[r, Mod[k, 4] + 1]] -> vid[[r, k]]}, {r, nR}, {k, 4}]], Union];
vdeg[v_] := Length[vnb[v]];
cwRays[v_] := SortBy[vnb[v], -ArcTan @@ (vpos[#] - vpos[v]) &];
stateTailId[s_] := vid[[ptStateRhomb[s], ptStateEdge[s] + 1]];

(* ---------- exit offsets at a vertex ---------- *)
Module[{v, cands, rays, p, R = 1.9, loc, d, items, names, a0},
  cands = Select[Keys[vnb], vdeg[#] == 7 && Norm[vpos[#]] < 12 &];
  v = First[SortBy[cands, Norm[vpos[#]] &]]; p = vpos[v]; rays = cwRays[v]; d = Length[rays];
  loc = near[p, R + 1.8];
  names = {"0 (戻る光線)\n(S,R) (L,R) (R,L)", "1\n(R,S) (S,L) (L,L)", "2\n(S,S) (L,S)", "3", "4", "5", "−1 ≡ 6\n(R,R)"};
  a0 = ArcTan @@ (vpos[rays[[1]]] - p);
  items = Table[With[{q = vpos[rays[[k]]]}, {Which[k == 1, Directive[GrayLevel[0.2], AbsoluteThickness[3.2]], k == 7, Directive[RGBColor[0.55, 0.3, 0.75], AbsoluteThickness[2.4]], k <= 3, Directive[ptRed, AbsoluteThickness[2.4]], True, Directive[GrayLevel[0.55], AbsoluteThickness[1.6]]], CapForm["Round"],
      If[k == 1, Arrow[{q, p + 0.14 (q - p)}], Arrow[{p, p + 0.92 (q - p)}]],
      Black, tx[names[[k]], p + 1.5 (q - p), 9, LineSpacing -> {1.1, 0}]}], {k, d}];
  ptSave["exit-offsets", frameW[{bgTiles[loc], Arrowheads[0.04], items, {Black, PointSize[0.022], Point[p]},
      {GrayLevel[0.3], AbsoluteThickness[1.], Arrowheads[0.03], Arrow[Table[p + 0.45 {Cos[a], Sin[a]}, {a, a0 - 0.3, a0 - 2.0, -0.1}]]},
      tx["時計回り", p + 0.72 {Cos[a0 - 1.15], Sin[a0 - 1.15]}, 8, GrayLevel[0.2]]}, p, 2.6, 2.0, ImageSize -> 430]]];

(* ---------- orbit search helpers on this window ---------- *)
depth = With[{x = cent[[All, 1]], y = cent[[All, 2]]}, MapThread[Min, {x - Min[x], Max[x] - x, y - Min[y], Max[y] - y}]];
classesOf[prog_String, d_ : 6.] := Module[{starts = Flatten[Table[st[r, e], {r, Flatten[Position[depth, _?(# > d &)]]}, {e, 0, 3}]], res},
   res = ptScanC[nav["sNext"], ptProg[prog], starts, 20000];
   ptClasses[nav, cent, ptProg[prog], starts, res]];
pivotPath[prog_String, s0_] := Module[{run = ptRun[nav["sNext"], ptProg[prog], s0]}, Append[First[edgeXY[#]] & /@ run["states"], First[edgeXY[s0]]]];
(* net walk: drop the excursion of every syntactic bounce (S followed by R) *)
netWalk[prog_String, s0_] := Module[{p = ptProg[prog], run, n, tails, keep, np},
   np = Length[p]; run = ptRun[nav["sNext"], p, s0]; n = Length[run["states"]];
   tails = stateTailId /@ run["states"];
   keep = Table[! (p[[Mod[k - 1, np] + 1]] == 1 && p[[Mod[k - 2, np] + 1]] == 0) && ! (p[[Mod[k - 1, np] + 1]] == 0 && p[[Mod[k, np] + 1]] == 1), {k, n}];
   Pick[RotateLeft[tails], keep]];
clsSR = classesOf["SR"]; clsSSR = classesOf["SSR"]; clsSSSR = classesOf["SSSR"]; cls6 = classesOf["SRSRSRSSR"];
pick[cls_, per_, nc_] := SelectFirst[cls, #period == per && #ncells == nc &];

(* ---------- pivot walks ---------- *)
Module[{panel},
  panel[prog_, c_, ttl_] := ptPanel[rh, xy, cent, c["cells"], "Margin" -> 0.9,
    "Extra" -> {{ptRed, AbsoluteThickness[1.3], JoinForm["Round"], Line[pivotPath[prog, c["rep"]]]}, {Black, PointSize[0.014], Point[DeleteDuplicates[pivotPath[prog, c["rep"]]]]}}, "Label" -> ttl];
  ptSave["pivot-walks", GraphicsRow[{
     panel["SR", pick[clsSR, 7, 7], "(SR)*: 周期 7"], panel["SSR", pick[clsSSR, 10, 10], "(SSR)*: 周期 10"],
     panel["SSR", pick[clsSSR, 28, 24], "(SSR)*: 周期 28"], panel["SSR", pick[clsSSR, 100, 70], "(SSR)*: 周期 100"]}, Spacings -> 4, ImageSize -> 600]]];

(* ---------- net exit sequences ---------- *)
Module[{panel},
  panel[prog_, c_, ttl_] := Module[{nw = netWalk[prog, c["rep"]], pts},
    pts = vpos /@ nw;
    ptPanel[rh, xy, cent, c["cells"], "Margin" -> 0.9, "Color" -> RGBColor[0.72, 0.84, 0.98],
     "Extra" -> {{ptGreen, AbsoluteThickness[1.4], JoinForm["Round"], Arrowheads[{{0.04, 0.55}}], Arrow /@ Partition[Append[pts, First[pts]], 2, 1]}, {Black, PointSize[0.018], Point[DeleteDuplicates[pts]]}}, "Label" -> ttl]];
  ptSave["net-walks", GraphicsRow[{
     panel["SSR", pick[clsSSR, 12, 14], "(3)* = (SSR)*: 周期 12"], panel["SSR", pick[clsSSR, 30, 35], "(3)* = (SSR)*: 周期 30"],
     panel["SSSR", clsSSSR[[2]], "(3,2)* = (SSSR)*: 周期 " <> ToString[clsSSSR[[2]]["period"]]],
     panel["SRSRSRSSR", Last[cls6], "(6)* = (SR)³(SSR): 周期 " <> ToString[Last[cls6]["period"]]]}, Spacings -> 4, ImageSize -> 600]]];

(* ---------- U-turn families ---------- *)
Module[{s = st[c0, 1], run1, run2, ctr, loc, arrows, p1 = "SSSRR", p2 = "SSSRRSRR", rb, panel},
  rb = ribbon[s, 9, 2]; ctr = cent[[ptStateRhomb[rb[[6]]]]]; loc = near[ctr, 8];
  run1 = ptRun[nav["sNext"], ptProg[p1], s, 2];
  run2 = ptRun[nav["sNext"], ptProg[p2], s, 3];
  arrows[run_] := Module[{pts = cent[[run["rhombs"]]], segs},
    segs = Partition[pts, 2, 1];
    Map[With[{dv = #[[2]] - #[[1]]}, With[{nn = 0.16 {-dv[[2]], dv[[1]]}/Norm[dv]}, {If[dv[[1]] > 0, ptBlue, ptRed], Arrow[{#[[1]] + nn + 0.12 dv, #[[2]] + nn - 0.12 dv}]}]] &, segs]];
  panel[run_, ttl_] := frameW[{bgTiles[loc], {RGBColor[1., 0.9, 0.7], EdgeForm[{GrayLevel[0.45], AbsoluteThickness[0.4]}], Polygon[xy[[DeleteDuplicates[run["rhombs"]]]]]},
     {AbsoluteThickness[1.4], Arrowheads[0.035], arrows[run]}}, ctr + {0.6, 0}, 4.4, 2.1, PlotLabel -> title[ttl]];
  ptSave["uturn", GraphicsRow[{panel[run1, "(S³RR)*: 周期 2 (往復)"], panel[run2, "(S³RR S RR)*: 1 サイクルごとに 2 タイル前進"]}, Spacings -> 8, ImageSize -> 600]]];
