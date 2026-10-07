(* corona_margin_exact.wl : exact value of delta_C = distance from the centre triangle t_C to the boundary of the
   union of the corona C, for the 30 coronas of the atlas, computed in exact algebraic arithmetic (no floating point):
   vertices are elements of Z[zeta] embedded with exact Cos/Sin of multiples of 72 degrees, boundary edges are
   found by exact matching of vertex keys, and squared distances are compared exactly (RootReduce). *)
keys = Import[FileNameJoin[{$dataDir, "atlas.wxf"}]];
e4x = Table[{Cos[2 Pi k/5], Sin[2 Pi k/5]}, {k, 0, 3}];                       (* exact basis 1, zeta, zeta^2, zeta^3 *)
embX[v_] := v . e4x;
sgn[x_] := With[{r = RootReduce[x]}, If[r === 0, 0, Sign[N[r, 60]]]];
segDist2[p_, {a_, b_}] := Module[{d = b - a, t, tt},
   t = ((p - a) . d)/(d . d);
   Which[sgn[t] <= 0, (p - a) . (p - a), sgn[t - 1] >= 0, (p - b) . (p - b), True, (p - a - t d) . (p - a - t d)]];
res = Table[Module[{t = ptCanonToTri[k], c = ptCanonCentre[k], tris, edges, bnd, cv, cands},
    tris = Table[{t["A"][[j]], t["B"][[j]], t["C"][[j]]}, {j, ptNTri[t]}];                    (* exact 4-vectors *)
    edges = Flatten[Table[Sort /@ Partition[Append[tris[[j]], tris[[j, 1]]], 2, 1], {j, Length[tris]}], 1];
    bnd = Select[Tally[edges], #[[2]] == 1 &][[All, 1]];                                       (* boundary edges, exact *)
    cv = embX /@ tris[[c]];
    cands = Join[Flatten[Table[segDist2[p, embX /@ e], {p, cv}, {e, bnd}]],
       Flatten[Table[segDist2[embX[q], {cv[[i]], cv[[Mod[i, 3] + 1]]}], {q, Union[Flatten[bnd, 1]]}, {i, 3}]]];
    First[SortBy[cands, N[RootReduce[#], 40] &]]], {k, keys}];
mins = RootReduce /@ res;
Print["delta_C^2 (exact) for the 30 coronas: ", Union[mins]];
Print["minimum delta^2 = ", RootReduce[First[SortBy[mins, N[#, 40] &]]], "  equals sin^2 36 = (5 - sqrt 5)/8: ", RootReduce[First[SortBy[mins, N[#, 40] &]] - (5 - Sqrt[5])/8] === 0];
Print["all delta_C >= 0.58 (exact): ", And @@ (sgn[# - (29/50)^2] >= 0 & /@ mins), "   numerically delta = ", N[Sqrt[First[SortBy[mins, N[#, 40] &]]], 12]];
