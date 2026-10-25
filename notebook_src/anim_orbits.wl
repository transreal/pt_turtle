(* anim_orbits.wl : animations of a turtle sweeping its orbit until the orbit closes.
   Every tile is coloured by the number of times the turtle has entered it so far (one colour scale common
   to all animations); the panel on the right shows the program, the counters and the histogram of tiles by
   number of visits.  One frame per command at the start, then several commands per frame so that a video
   lasts about half a minute.  A spec may instead show one orbit of every congruence class in turn
   ("classes" -> True), and a large orbit may be followed by a zooming camera ("camera" -> "follow") with a
   geometric speed-up ("ramp" -> True).  "escape" -> True follows a turtle that never closes (diffusive regime)
   from the centre of the window until it reaches the border, with a log-log plot of its distance from the start.
   usage : wolfram -noinit -noprompt -script run.wl anim_orbits.wl [name ...]        (default: all)
           wolfram -noinit -noprompt -script run.wl anim_orbits.wl test [name ...]   (4 still frames per orbit)
   output: anim/turtle-<name>.mp4 (1920x1080, 30 fps, H.264) and anim/turtle-<name>-final.png.
   needs ffmpeg on the PATH or in the environment variable FFMPEG.  A window of side n is read from
   data/win<n>.wxf when present (bigwin.wl n), otherwise generated (side 300: about a minute). *)

animSpecs = {
   <|"name" -> "sr", "prog" -> "SR", "classes" -> True, "title" -> "(SR)*", "regime" -> "有界", "segLabel" -> "頂点のまわりの並び方",
     "slowAll" -> True, "startHold" -> 20, "segEndHold" -> 45, "minSpan" -> 7,
     "note" -> "S と R を交互に繰り返すと、タートルは 1 つの頂点のまわりを
一周して元に戻る。頂点に集まるタイルの枚数 (3〜7) がそのまま
周期で、どのタイルも 1 回だけ通る。並び方 7 通りを順に見せる。"|>,
   <|"name" -> "ssr", "prog" -> "SSR", "period" -> 100, "title" -> "(SSR)*", "regime" -> "有界",
     "note" -> "SR8 の最短ワームの頭の動き。どのペンローズ菱形タイリングでも
軌道は必ず閉じ、周期は 10, 12, 28, 30, 50, 80, 100 のどれか。
ここは最大の周期 100 の軌道。"|>,
   <|"name" -> "sr3ssr", "prog" -> "SRSRSRSSR", "period" -> 50, "title" -> "(SR)³(SSR)", "regime" -> "有界",
     "note" -> "有界なプログラムの一つ。S と R が交互に続く間は狭い範囲を
回り続けるので、同じタイルに何度も入る (最大 18 回)。
周期 50 の軌道。"|>,
   <|"name" -> "ssr-ssssr", "prog" -> "SSRSSSSR", "period" -> 309, "title" -> "(SSR)(SSSSR)", "regime" -> "有界",
     "note" -> "有界。軌道は合同なものを同一視して 13 種類あり、
ここは最大の周期 309 の軌道。中心に近いタイルほど何度も通る。"|>,
   <|"name" -> "sssr", "prog" -> "SSSR", "period" -> 1385, "title" -> "(SSSR)*", "regime" -> "非有界 (階層型)",
     "note" -> "非有界。同じ中心をもつ環状の軌道が入れ子に並び、外側ほど
周期が長い (周期に上限がない)。ここは周期 1385 の環状の軌道。"|>,
   <|"name" -> "sslsr", "prog" -> "SSLSR", "period" -> 370, "title" -> "(SSLSR)*", "regime" -> "拡散型",
     "note" -> "ほとんどの軌道は閉じずに遠くへ広がっていく。一辺 90 の窓で
調べた範囲で閉じる軌道は、この周期 370 の類だけ。"|>,
   <|"name" -> "sr2sssr", "prog" -> "SRSRSSSR", "period" -> 16265, "window" -> 300, "camera" -> "follow", "ramp" -> True,
     "title" -> "(SR)²(SSSR)", "regime" -> "非有界 (階層型)",
     "note" -> "雪の結晶 (コッホ曲線) のように何段にも折れ曲がった軌道。
窓を広げるとさらに大きな軌道が現れ、周期に上限はない。一辺 300 の
窓で最大の周期 16265 の軌道。カメラはタートルを追って引いていく。"|>,
   <|"name" -> "sslsr-walk", "prog" -> "SSLSR", "escape" -> True, "window" -> 300, "camera" -> "follow", "ramp" -> True,
     "title" -> "(SSLSR)*", "regime" -> "拡散型",
     "note" -> "閉じない軌道の典型例。窓の中心から出発し、行きつ戻りつしながら
ランダムウォークのように広がって一辺 300 の窓の端に達する。出発点
からの距離はおおよそ時間の平方根に比例して増える (下のグラフ)。"|>};

fps = 30; imgW = 1920; imgH = 1080; side = 1080;
$animDir = FileNameJoin[{$buildDir, "anim"}];
If[! DirectoryQ[$animDir], CreateDirectory[$animDir]];

ffmpegExe = Module[{e = Environment["FFMPEG"], w},
   If[StringQ[e] && FileExistsQ[e], e,
    w = Quiet[RunProcess[{If[$OperatingSystem === "Windows", "where", "which"], "ffmpeg"}, "StandardOutput"]];
    If[StringQ[w] && StringTrim[w] =!= "", StringTrim[First[StringSplit[StringTrim[w], "\n"]]], $Failed]]];
Print["ffmpeg: ", ffmpegExe];

(* ---- pentagrid windows (cached per session); sets rh, nav, sNext, xy, cent, fat, nTiles, starts, ctr ---- *)
$winCache = <||>;
useWindow[n_Integer] := Module[{f = FileNameJoin[{$dataDir, "win" <> ToString[n] <> ".wxf"}], t0 = AbsoluteTime[], w, r, sn, xyw, ce, inter},
   If[! KeyExistsQ[$winCache, n],
    If[FileExistsQ[f],
     w = Import[f]; r = <|"quad" -> w["quad"], "fat" -> w["fat"]|>; sn = w["sNext"]; xyw = w["xy"]; ce = w["cent"],
     r = ptPentagrid[{-3, -3}, {n, n}]; sn = ptNav[r]["sNext"]; xyw = Developer`ToPackedArray[r["quad"] . ptE4]; ce = Total[xyw, {2}]/4.];
    inter = Flatten[Position[Min /@ Partition[sn, 4], _?Positive]];
    $winCache[n] = <|"rh" -> r, "sNext" -> sn, "xy" -> xyw, "cent" -> ce, "starts" -> Flatten[Outer[Plus, 4 (inter - 1), {1, 2, 3, 4}]]|>;
    Print["window ", n, ": ", Length[xyw], " tiles  [", Round[AbsoluteTime[] - t0], "s]"]];
   w = $winCache[n];
   rh = w["rh"]; sNext = w["sNext"]; nav = <|"sNext" -> sNext|>; xy = w["xy"]; cent = w["cent"]; fat = rh["fat"]; nTiles = Length[xy];
   starts = w["starts"]; ctr = Mean[cent]; n];

(* closed orbit of the given period nearest to the window centre; returns its start state (program position 0) *)
pickOrbit[prog_String, period_Integer] := Module[{p = ptProg[prog], res, oids},
   res = ptScanC[sNext, p, starts, Max[20000, period + 100]];
   oids = DeleteDuplicates[Select[res, #[[1]] == 0 && #[[2]] == period &][[All, 3]]];
   If[oids === {}, $Failed,
    First[SortBy[oids, Norm[Mean[cent[[DeleteDuplicates[ptOrbitRhombsC[sNext, p, #, period]]]]] - ctr] &]]]];

(* the orbits of a spec: {{start state, period}, ...} *)
orbitList[spec_] := Module[{p = ptProg[spec["prog"]], st, res, cls, s0, c0, runs},
   If[TrueQ[spec["escape"]],
    (* the tile nearest to the window centre; of its four headings, the one whose walk reaches the border last *)
    c0 = First[MinimalBy[DeleteDuplicates[ptStateRhomb /@ starts], Norm[cent[[#]] - ctr] &]];
    runs = Table[{4 (c0 - 1) + e + 1, ptRunC[sNext, p, 4 (c0 - 1) + e + 1, 10^6]}, {e, 0, 3}];
    Print[spec["name"], ": centre tile ", c0, ", {start, {flag (1 = escaped), cycles}} = ", runs];
    runs = Select[runs, #[[2, 1]] == 1 &];
    Return[If[runs === {}, $Failed, {{#[[1]], #[[2, 2]], "escape"} &[First[MaximalBy[runs, #[[2, 2]] &]]]}]]];
   If[TrueQ[spec["classes"]],
    (* one orbit of every congruence class, taken from the starts near the window centre, by period *)
    st = Select[starts, Norm[cent[[ptStateRhomb[#]]] - ctr] < 12 &];
    res = ptScanC[sNext, p, st, 20000];
    cls = ptExactClasses[rh["quad"], sNext, cent, p, res];
    Print[spec["name"], ": ", Length[cls], " classes {period, tiles}: ", {#period, #ncells} & /@ cls];
    {#rep, #period} & /@ cls,
    s0 = pickOrbit[spec["prog"], spec["period"]];
    If[s0 === $Failed, $Failed, {{s0, spec["period"]}}]]];

(* states after t = 0 .. T commands, and the tile entered by command t (0 for a turn) *)
recordOrbit[prog_String, s0_Integer, period_Integer] := Module[{p = ptProg[prog], np, T, st, ent, s = s0, c},
   np = Length[p]; T = np period;
   st = ConstantArray[0, T + 1]; ent = ConstantArray[0, T]; st[[1]] = s0;
   Do[c = p[[Mod[t - 1, np] + 1]]; s = ptApply[nav, s, c]; st[[t + 1]] = s; If[c == 0, ent[[t]] = ptStateRhomb[s]], {t, T}];
   If[s != s0, Print["  orbit does not close after ", period, " cycles!"]];
   <|"prog" -> p, "np" -> np, "T" -> T, "states" -> st, "enter" -> ent|>];

(* a walk that leaves the window: the commands up to (not including) the first S move without a neighbour *)
recordWalk[prog_String, s0_Integer, maxCyc_Integer] := Module[{p = ptProg[prog], np, st, ent, s = s0, c, t = 0, s2},
   np = Length[p]; st = Internal`Bag[{s0}]; ent = Internal`Bag[];
   Catch[Do[c = p[[j]];
      If[c == 0, s2 = sNext[[s]]; If[s2 == 0, Throw[Null]]; s = s2, s = ptApply[nav, s, c]];
      t++; Internal`StuffBag[st, s]; Internal`StuffBag[ent, If[c == 0, ptStateRhomb[s], 0]], {cyc, maxCyc + 1}, {j, np}]];
   <|"prog" -> p, "np" -> np, "T" -> t, "states" -> Internal`BagPart[st, All], "enter" -> Internal`BagPart[ent, All]|>];

(* rendered states: {t, number of video frames showing the state after t commands, commands per frame} *)
schedule[T_, np_, spec_, last_] := Module[{hold = Lookup[spec, "hold", 6], p1, p2, pts, budget = 22 fps, r = 1.06, k, kMax, t, j, framesFor, lo, hi},
   p1 = If[TrueQ[spec["slowAll"]], T, Min[T, 3 np]]; p2 = Min[T, p1 + 6 np];
   pts = Join[{{0, Lookup[spec, "startHold", fps], 0}}, Table[{t, hold, 1}, {t, 1, p1}], Table[{t, 2, 1}, {t, p1 + 1, p2}]];
   If[p2 < T,
    If[TrueQ[spec["ramp"]],
     (* commands per frame grow geometrically up to kMax, the smallest value that fits the rest into the budget *)
     framesFor[km_] := Module[{tt = p2, jj = 0, f = 0}, While[tt < T, jj++; tt += Min[km, Ceiling[r^jj]]; f++]; f];
     lo = 1; hi = T - p2;
     While[lo < hi, With[{mid = Quotient[lo + hi, 2]}, If[framesFor[mid] <= budget, hi = mid, lo = mid + 1]]];
     kMax = lo; t = p2; j = 0;
     While[t < T, j++; k = Min[kMax, Ceiling[r^j]]; t = Min[T, t + k]; AppendTo[pts, {t, 1, k}]],
     k = Max[1, Ceiling[(T - p2)/budget]];
     pts = Join[pts, Table[{t, 1, k}, {t, p2 + k, T, k}]];
     If[pts[[-1, 1]] != T, AppendTo[pts, {T, 1, k}]]]];
   pts[[-1, 2]] = If[last, 4 fps, Lookup[spec, "segEndHold", 2 fps]];
   pts];

(* square view {cx, cy, span} around a bounding box *)
boxView[{x0_, x1_}, {y0_, y1_}, minSpan_ : 0] := With[{s = Max[x1 - x0, y1 - y0]}, {(x0 + x1)/2, (y0 + y1)/2, Max[minSpan, s + 2 (1.4 + 0.03 s)]}];
(* shift a view so that the point c stays inside *)
clampView[{cx_, cy_, span_}, c_] := Module[{h = span/2 - 1.5, x = cx, y = cy},
   If[c[[1]] > x + h, x = c[[1]] - h]; If[c[[1]] < x - h, x = c[[1]] + h];
   If[c[[2]] > y + h, y = c[[2]] - h]; If[c[[2]] < y - h, y = c[[2]] + h];
   {x, y, span}];

(* colour of a tile entered k times: yellow, orange, red, crimson, magenta, purple, indigo, blue, ... darkening to navy for k >= 20 *)
visitCols = RGBColor /@ {"#FFE066", "#FBA43C", "#EB5A2A", "#C81D4E", "#A0187A", "#74208F", "#4B2A99", "#2E3D98", "#22508A", "#1B3D6B", "#142C4F"};
visitCol[k_] := If[k <= Length[visitCols], visitCols[[k]], Blend[{Last[visitCols], RGBColor["#070D1A"]}, Min[1., (k - Length[visitCols])/9.]]];
turtleRed = RGBColor[0.86, 0.08, 0.10]; trailCol = RGBColor["#0B8FA6"];
jp[s_, size_, opts___] := Style[s, FontFamily -> ptFont, FontSize -> size, opts];
num[s_, size_, opts___] := Style[s, FontFamily -> "Consolas", FontSize -> size, opts];

(* frames of one orbit, appended to frameDir from number fno0 + 1 on; returns the last frame number *)
renderOrbit[spec_, {s0_, period_, mode_ : "closed"}, seg_, nSeg_, frameDir_, fno0_, testOnly_, tStart_] := Module[{
    escape = (mode === "escape"), distS, dMax, aRef, lay, distPlot,
    name = spec["name"], rec, p, np, T, st, ent, tile0, lastEnt, entIdx, entTiles, nS, cntF, K, nCells, sched, holdOf, kOf, idxOf,
    last = (seg == nSeg), follow = (Lookup[spec, "camera", None] === "follow"), cellsAll, pts, view0, cand, views, bb, rawV, perFrame, smooth,
    cnt, hist, hMax, m, fno = fno0, tileG, panelG, frameImg, img, hx0 = 1140., hw = 720., hb = 190., ht = 285., tests},
   rec = If[escape, recordWalk[spec["prog"], s0, period], recordOrbit[spec["prog"], s0, period]];
   p = rec["prog"]; np = rec["np"]; T = rec["T"]; st = rec["states"]; ent = rec["enter"];
   tile0 = ptStateRhomb[s0];
   entIdx = Flatten[Position[ent, _?Positive]]; entTiles = ent[[entIdx]]; lastEnt = If[escape, 0, Last[entIdx]]; nS = Length[entIdx];
   (* final visit counts: the last S move re-enters the start tile and closes the first stay there *)
   cntF = ConstantArray[0, nTiles]; cntF[[tile0]] = 1;
   Do[If[t != lastEnt, cntF[[ent[[t]]]] += 1], {t, entIdx}];
   K = Max[cntF]; nCells = Count[cntF, _?Positive];
   sched = schedule[T, np, spec, last];
   holdOf = Association[#[[1]] -> #[[2]] & /@ sched]; kOf = Association[#[[1]] -> #[[3]] & /@ sched];
   idxOf = AssociationThread[sched[[All, 1]] -> Range[Length[sched]]];
   Print[name, If[nSeg > 1, " [" <> ToString[seg] <> "/" <> ToString[nSeg] <> "]", ""], ": (", spec["prog"], ")* ", If[escape, "reaches the border after " <> ToString[Quotient[T, np]] <> " cycles", "period " <> ToString[period]], ", ", T, " commands, ", nS, " S moves, ", nCells,
    " tiles, max visits ", K, ", visit tally ", Sort[Tally[Select[cntF, Positive]]], "; ", Length[sched], " rendered frames, ", Total[sched[[All, 2]]], " video frames, at most ", Max[sched[[All, 3]]], " commands per frame"];

   (* the view of the whole orbit; background candidates around it *)
   cellsAll = Flatten[Position[cntF, _?Positive]];
   pts = Flatten[xy[[cellsAll]], 1];
   view0 = boxView[MinMax[pts[[All, 1]]], MinMax[pts[[All, 2]]], Lookup[spec, "minSpan", 0]];
   cand = With[{h = view0[[3]]/2 + If[follow, 20, 2]},
     Flatten[Position[UnitStep[cent[[All, 1]] - view0[[1]] + h] UnitStep[view0[[1]] + h - cent[[All, 1]]] UnitStep[cent[[All, 2]] - view0[[2]] + h] UnitStep[view0[[2]] + h - cent[[All, 2]]], 1]]];

   tileG[t_, groups_, trail_, view_] := Module[{s = st[[t + 1]], r, de, c, mm, x0, x1, y0, y1, ppu, edgeTh, curTh, headTh, arTh, ah, ringR, bg, thinCol, bgEdge, paintEdge},
     r = ptStateRhomb[s]; de = ptStateDExy[rh, s]; c = cent[[r]]; mm = Mean[de];
     {x0, x1} = view[[1]] + {-1, 1} view[[3]]/2; {y0, y1} = view[[2]] + {-1, 1} view[[3]]/2;
     ppu = side/view[[3]];                              (* pixels per unit length *)
     edgeTh = Clip[ppu/110, {0.25, 0.8}]; curTh = Clip[ppu/22, {1.2, 3.5}]; headTh = Clip[ppu/11, {2.2, 7.}]; arTh = Clip[ppu/40, {0.7, 2.6}];
     ah = Clip[0.28 ppu/side, {0.004, 0.03}]; ringR = If[ppu < 40, 24./ppu, None];
     (* level of detail for the zooming camera: thin tiles fade into the background, edges fade out *)
     bg = If[follow && ppu < 8, {},
       Pick[cand, UnitStep[cent[[cand, 1]] - x0 + 1.5] UnitStep[x1 + 1.5 - cent[[cand, 1]]] UnitStep[cent[[cand, 2]] - y0 + 1.5] UnitStep[y1 + 1.5 - cent[[cand, 2]]], 1]];
     thinCol = If[follow, Blend[{ptColFat, ptColThin}, Clip[(ppu - 8)/4, {0, 1}]], ptColThin];
     bgEdge = Which[! follow || ppu >= 20, EdgeForm[{GrayLevel[0.76], AbsoluteThickness[edgeTh]}], ppu >= 12, EdgeForm[{GrayLevel[0.76 + 0.2 (20 - ppu)/8], AbsoluteThickness[edgeTh]}], True, EdgeForm[None]];
     paintEdge = Which[! follow || ppu >= 20, EdgeForm[{GrayLevel[0.38], AbsoluteThickness[edgeTh]}], ppu >= 10, EdgeForm[{GrayLevel[0.38], Opacity[(ppu - 10)/10], AbsoluteThickness[edgeTh]}], True, EdgeForm[None]];
     Graphics[{
       If[follow, {ptColFat, Rectangle[{x0, y0}, {x1, y1}]}, {}],
       If[bg === {}, {}, {bgEdge, {ptColFat, Polygon[xy[[Pick[bg, fat[[bg]], 1]]]]}, {thinCol, Polygon[xy[[Pick[bg, fat[[bg]], 0]]]]}}],
       {paintEdge, {visitCol[#[[1]]], Polygon[xy[[#[[2]]]]]} & /@ groups},
       {{White, AbsoluteThickness[Max[1.2, curTh/2] + 2], Circle[cent[[tile0]], Max[0.2, 7./ppu]]}, {Black, AbsoluteThickness[Max[1.2, curTh/2]], Circle[cent[[tile0]], Max[0.2, 7./ppu]]}},
       If[Length[trail] > 1, {CapForm["Round"], JoinForm["Round"], {White, Opacity[0.85], AbsoluteThickness[Clip[ppu/30, {1.2, 2.5}] + 2.5], Line[cent[[trail]]]}, {trailCol, AbsoluteThickness[Clip[ppu/30, {1.2, 2.5}]], Line[cent[[trail]]]}}, {}],
       {FaceForm[], EdgeForm[{Black, AbsoluteThickness[curTh]}], Polygon[xy[[r]]]},
       {turtleRed, AbsoluteThickness[headTh], CapForm["Round"], Line[de]},
       {Black, AbsoluteThickness[arTh], Arrowheads[ah], Arrow[{c - 0.15 (mm - c), c + 0.72 (mm - c)}]},
       If[ringR =!= None, {{White, AbsoluteThickness[5.5], Circle[c, ringR]}, {turtleRed, AbsoluteThickness[2.5], Circle[c, ringR]}}, {}]},
      PlotRange -> {{x0, x1}, {y0, y1}}, PlotRangePadding -> 0, ImagePadding -> 0, AspectRatio -> 1, PlotRangeClipping -> True, Background -> White]];

   lay = If[escape, <|"hb" -> 408., "ht" -> 70., "lab" -> 72|>, <|"hb" -> hb, "ht" -> ht, "lab" -> 80|>];
   panelG[t_, hist_, nPainted_, nMoves_] := Module[{j = If[t == 0, 0, Mod[t - 1, np] + 1], b = 52., gap = 8., slot, gx0, done = (t == T), v, hb = lay["hb"], ht = lay["ht"], tot, ys = If[escape, 40, 48]},
     slot = Min[hw/K, 110.]; gx0 = hx0 + (hw - slot K)/2;
     tot[x_] := If[escape, "", " / " <> ToString[x]];
     {
      Text[jp[spec["title"], 54, Bold, Black], {hx0, 1018}, {-1, 0}],
      Text[jp[spec["regime"] <> If[nSeg > 1, "      " <> Lookup[spec, "segLabel", "軌道"] <> " " <> ToString[seg] <> " / " <> ToString[nSeg], ""], 26, Bold, GrayLevel[0.35]], {hx0, 962}, {-1, 0}],
      Inset[Pane[jp[spec["note"], 20, GrayLevel[0.25], LineSpacing -> {1.25, 0}, TextAlignment -> Left], {hw, Automatic}], {hx0, 935}, {Left, Top}],
      Text[jp["命令列 (先頭に戻って繰り返す)", 20, GrayLevel[0.35]], {hx0, 797}, {-1, 0}],
      Table[With[{xl = hx0 + (i - 1) (b + gap), on = (i == j)},
        {EdgeForm[{GrayLevel[0.45], AbsoluteThickness[1]}], If[on, turtleRed, White], Rectangle[{xl, 727}, {xl + b, 727 + b}],
         Text[jp[{"S", "R", "L"}[[p[[i]] + 1]], 32, Bold, If[on, White, Black]], {xl + b/2, 727 + b/2}]}], {i, np}],
      Text[Row[{jp["サイクル  ", 28, GrayLevel[0.3]], num[ToString[Quotient[t, np]] <> tot[period], 30, Black]}], {hx0, 665}, {-1, 0}],
      Text[Row[{jp["S 移動  ", 28, GrayLevel[0.3]], num[ToString[nMoves] <> tot[nS], 30, Black]}], {hx0, 622}, {-1, 0}],
      Text[Row[{jp["塗ったタイル  ", 28, GrayLevel[0.3]], num[ToString[nPainted] <> tot[nCells], 30, Black]}], {hx0, 579}, {-1, 0}],
      Text[jp["タイルの色 = これまでにそのタイルに入った回数", 22, Black], {hx0, 528}, {-1, 0}],
      {GrayLevel[0.4], AbsoluteThickness[1], Line[{{gx0 - 6, hb}, {gx0 + slot K + 6, hb}}]},
      Table[With[{xl = gx0 + (k - 1) slot + 0.12 slot, xr = gx0 + k slot - 0.12 slot, h = ht hist[[k]]/hMax},
        {If[hist[[k]] > 0, {EdgeForm[{GrayLevel[0.3], AbsoluteThickness[0.8]}], visitCol[k], Rectangle[{xl, hb}, {xr, hb + Max[h, 1.5]}],
           Text[num[ToString[hist[[k]]], If[K > 12, 15, 18], GrayLevel[0.2]], {(xl + xr)/2, hb + Max[h, 1.5] + 14}]}, {}],
         {EdgeForm[{GrayLevel[0.6], AbsoluteThickness[0.6]}], visitCol[k], Rectangle[{xl, hb - 30}, {xr, hb - 10}]},
         Text[num[ToString[k], If[K > 12, 15, 19], Black], {(xl + xr)/2, hb - 46}]}], {k, K}],
      Text[jp["通過回数", 20, GrayLevel[0.3]], {hx0 + hw/2, hb - lay["lab"]}],
      If[escape, distPlot[t], {}],
      Which[
       t == 0, Text[jp["黒枠 = タートルのいるタイル、赤い辺 = 進行方向、○ = 出発タイル、青緑の線 = 直前の足跡", 19, GrayLevel[0.25]], {hx0, ys}, {-1, 0}],
       done && escape, Text[jp["窓の端に達した: " <> ToString[Quotient[T, np]] <> " サイクルたっても閉じずに広がり続けた", 25, Bold, turtleRed], {hx0, ys}, {-1, 0}],
       done, Text[jp["元の状態に戻った: 周期 " <> ToString[period] <> " サイクル (命令 " <> ToString[T] <> " 回)", 27, Bold, turtleRed], {hx0, ys}, {-1, 0}],
       True, v = fps kOf[t]/holdOf[t];
       Text[jp["速さ: 毎秒 " <> ToString[If[IntegerQ[v], v, Round[N[v], 0.1]]] <> " 命令", 21, GrayLevel[0.45]], {hx0, ys}, {-1, 0}]]}];

   (* escape mode: distance of the turtle from the start tile against cycles, log-log, with a reference line of slope 1/2 *)
   distPlot[t_] := Module[{x0p = hx0 + 64, x1p = hx0 + hw - 24, y0p = 112., y1p = 262., lc, ld, fx, fy, cur = distS[[1 ;; idxOf[t]]], pl, cM = Max[10., T/np]},
     lc = Log10[cM]; ld = Log10[Max[10., 1.15 Max[dMax, aRef Sqrt[cM]]]];
     fx[c_] := x0p + (x1p - x0p) Log10[Max[1., c]]/lc; fy[d_] := y0p + (y1p - y0p) Log10[Max[1., d]]/ld;
     pl = {fx[#[[1]]], fy[#[[2]]]} & /@ Select[cur, #[[1]] >= 1 &];
     {Text[jp["出発点からの距離 (両対数目盛)", 22, Black], {hx0, 294}, {-1, 0}],
      {GrayLevel[0.4], AbsoluteThickness[1], Line[{{x0p, y1p}, {x0p, y0p}, {x1p, y0p}}]},
      Table[{{GrayLevel[0.4], Line[{{fx[10^i], y0p}, {fx[10^i], y0p - 6}}]}, Text[num[{"1", "10", "100", "1000", "10⁴", "10⁵", "10⁶"}[[i + 1]], 15, GrayLevel[0.3]], {fx[10^i], y0p - 20}]}, {i, 0, Min[6, Floor[lc]]}],
      Table[{{GrayLevel[0.4], Line[{{x0p, fy[10^i]}, {x0p - 6, fy[10^i]}}]}, Text[num[{"1", "10", "100", "1000"}[[i + 1]], 15, GrayLevel[0.3]], {x0p - 12, fy[10^i]}, {1, 0}]}, {i, 0, Min[3, Floor[ld]]}],
      Text[jp["サイクル数", 16, GrayLevel[0.35]], {x1p, y0p - 42}, {1, 0}],
      {GrayLevel[0.55], Dashing[{6, 5}], AbsoluteThickness[1.4], Line[{{fx[1.], fy[aRef]}, {fx[cM], fy[aRef Sqrt[cM]]}}]},
      Text[jp["破線: 傾き 1/2", 16, GrayLevel[0.4]], {x0p + 12, y1p - 4}, {-1, 1}],
      If[Length[pl] > 1, {trailCol, AbsoluteThickness[1.6], Line[pl]}, {}],
      If[pl =!= {}, {turtleRed, Disk[Last[pl], 5]}, {}]}];

   frameImg[t_, hist_, nPainted_, nMoves_, groups_, trail_, view_] := Module[{g, im},
     g = Graphics[{Inset[tileG[t, groups, trail, view], {side/2, side/2}, Center, {side, side}],
        {GrayLevel[0.82], AbsoluteThickness[1], Line[{{side, 0}, {side, imgH}}]},
        panelG[t, hist, nPainted, nMoves]},
       PlotRange -> {{0, imgW}, {0, imgH}}, ImageSize -> {imgW, imgH}, AspectRatio -> imgH/imgW, ImagePadding -> 0, PlotRangePadding -> 0, Background -> White];
     im = Rasterize[g, ImageResolution -> 72, Background -> White];
     If[ImageDimensions[im] =!= {imgW, imgH}, im = ImageResize[im, {imgW, imgH}]];
     im];

   (* walk through the orbit; pass 1 finds the scale of the histogram and the camera path, pass 2 renders *)
   tests = If[TrueQ[testOnly], Union[First[Nearest[sched[[All, 1]], #]] & /@ {0, Min[T, 2 np + 1], Quotient[T, 2], T}], None];
   Do[
    cnt = ConstantArray[0, nTiles]; cnt[[tile0]] = 1; hist = ConstantArray[0, K]; hist[[1]] = 1; m = 0;
    If[pass == 1, hMax = 1; rawV = {}; distS = {}; bb = Join[MinMax[xy[[tile0, All, 1]]], MinMax[xy[[tile0, All, 2]]]]];
    Do[
     If[t > 0 && ent[[t]] > 0,
      m++;
      If[t != lastEnt, With[{e = ent[[t]], c = cnt[[ent[[t]]]]}, If[c > 0, hist[[c]] -= 1]; cnt[[e]] = c + 1; hist[[c + 1]] += 1]];
      If[pass == 1 && follow, With[{q = xy[[ent[[t]]]]},
        bb = {Min[bb[[1]], Min[q[[All, 1]]]], Max[bb[[2]], Max[q[[All, 1]]]], Min[bb[[3]], Min[q[[All, 2]]]], Max[bb[[4]], Max[q[[All, 2]]]]}]]];
     If[pass == 1,
      hMax = Max[hMax, Max[hist]];
      If[follow && KeyExistsQ[holdOf, t], AppendTo[rawV, boxView[bb[[1 ;; 2]], bb[[3 ;; 4]], 14.]]];
      If[escape && KeyExistsQ[holdOf, t], AppendTo[distS, {N[t/np], Norm[cent[[ptStateRhomb[st[[t + 1]]]]] - cent[[tile0]]]}]],
      If[If[tests === None, KeyExistsQ[holdOf, t], MemberQ[tests, t]],
       Module[{nz = Pick[Range[nTiles], Unitize[cnt], 1], groups, trail, f, tl},
        groups = KeyValueMap[{#1, #2} &, GroupBy[nz, cnt[[#]] &]];
        tl = Max[20, Min[400, Ceiling[1.5 kOf[t] nS/T]]];
        trail = If[t == T, {}, Join[If[m < tl, {tile0}, {}], entTiles[[Max[1, m - tl + 1] ;; m]], {ptStateRhomb[st[[t + 1]]]}]];
        img = frameImg[t, hist, Length[nz], m, groups, trail, If[follow, clampView[views[[idxOf[t]]], cent[[ptStateRhomb[st[[t + 1]]]]]], view0]];
        If[tests =!= None,
         Export[FileNameJoin[{$animDir, "test-" <> name <> If[nSeg > 1, "-" <> ToString[seg], ""] <> "-" <> ToString[t] <> ".png"}], img];
         Print["  test frame t=", t, "  [", Round[AbsoluteTime[] - tStart, 0.1], "s]"],
         fno++; f = FileNameJoin[{frameDir, IntegerString[fno, 10, 6] <> ".png"}]; Export[f, img];
         Do[fno++; CopyFile[f, FileNameJoin[{frameDir, IntegerString[fno, 10, 6] <> ".png"}]], {holdOf[t] - 1}];
         If[t == T && last, Export[FileNameJoin[{$animDir, "turtle-" <> name <> "-final.png"}], img]];
         If[Mod[fno, 200] < holdOf[t], Print["  frame ", fno, " (t = ", t, " / ", T, ")  [", Round[AbsoluteTime[] - tStart], "s]"]]]]]],
     {t, 0, T}];
    If[pass == 1 && escape,
     dMax = Max[distS[[All, 2]]];
     aRef = With[{q = Select[distS, #[[1]] >= 100 && #[[2]] > 0 &]}, If[q === {}, 1., Exp[Mean[Log[q[[All, 2]]] - 0.5 Log[q[[All, 1]]]]]]]];
    (* camera: spans never shrink; centre and span are averaged over one second of video either side *)
    If[pass == 1 && follow,
     rawV[[All, 3]] = FoldList[Max, rawV[[All, 3]]];
     perFrame = Join @@ MapThread[ConstantArray[#1, #2] &, {rawV, sched[[All, 2]]}];
     smooth = Table[Mean[perFrame[[Max[1, i - fps] ;; Min[Length[perFrame], i + fps]]]], {i, Length[perFrame]}];
     views = smooth[[Most[Prepend[Accumulate[sched[[All, 2]]], 0]] + 1]]],
    {pass, 1, 2}];
   fno];

makeAnim[spec_, testOnly_ : False] := Module[{name = spec["name"], orbits, frameDir = None, fno = 0, out, tStart = AbsoluteTime[]},
   useWindow[Lookup[spec, "window", 90]];
   orbits = orbitList[spec];
   If[! ListQ[orbits] || orbits === {}, Print[name, ": no closed orbit found"]; Return[$Failed]];
   If[! TrueQ[testOnly],
    frameDir = FileNameJoin[{$TemporaryDirectory, "turtle-anim-" <> name}];
    If[DirectoryQ[frameDir], DeleteDirectory[frameDir, DeleteContents -> True]]; CreateDirectory[frameDir]];
   Do[fno = renderOrbit[spec, orbits[[i]], i, Length[orbits], frameDir, fno, testOnly, tStart], {i, Length[orbits]}];
   If[TrueQ[testOnly], Return[Null]];
   out = FileNameJoin[{$animDir, "turtle-" <> name <> ".mp4"}];
   If[ffmpegExe === $Failed, Print["  ffmpeg not found; frames left in ", frameDir]; Return[frameDir]];
   With[{r = RunProcess[{ffmpegExe, "-y", "-loglevel", "error", "-framerate", ToString[fps], "-start_number", "1", "-i", FileNameJoin[{frameDir, "%06d.png"}],
        "-c:v", "libx264", "-preset", "slow", "-crf", "20", "-pix_fmt", "yuv420p", "-movflags", "+faststart", out}]},
    If[r["ExitCode"] != 0, Print["  ffmpeg failed: ", r["StandardError"]]; Return[frameDir]]];
   DeleteDirectory[frameDir, DeleteContents -> True];
   Print["  wrote ", out, ": ", fno, " frames = ", Round[fno/fps, 0.1], " s, ", Round[FileByteCount[out]/2.^20, 0.1], " MB  [", Round[AbsoluteTime[] - tStart], "s]"];
   out];

$animArgs = Rest[$args];
$animTest = MemberQ[$animArgs, "test"];
$animNames = DeleteCases[$animArgs, "test"];
Do[If[$animNames === {} || MemberQ[$animNames, spec["name"]], makeAnim[spec, $animTest]], {spec, animSpecs}];
