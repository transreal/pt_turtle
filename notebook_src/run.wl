(* ASCII launcher: wolfram -noinit -noprompt -script run.wl target.wl [args...] ; loads UTF-8 sources *)
$CharacterEncoding = "UTF-8";
$buildDir = DirectoryName[$InputFileName];
(* repository root: this folder is <repo>/notebook_src *)
$repoDir = ParentDirectory[$buildDir] <> "/";
$figDir = FileNameJoin[{$buildDir, "figs"}];
$dataDir = FileNameJoin[{$buildDir, "data"}];
If[! DirectoryQ[$figDir], CreateDirectory[$figDir]];
If[! DirectoryQ[$dataDir], CreateDirectory[$dataDir]];
$args = With[{p = Position[$CommandLine, _String?(StringEndsQ[#, "run.wl"] &)]}, If[p === {}, {}, Drop[$CommandLine, p[[1, 1]]]]];
Get[FileNameJoin[{$buildDir, "lib.wl"}], CharacterEncoding -> "UTF-8"];
$t0 = AbsoluteTime[];
Get[FileNameJoin[{$buildDir, First[$args]}], CharacterEncoding -> "UTF-8"];
Print["[done ", First[$args], " ", Round[AbsoluteTime[] - $t0, 0.1], "s]"];
Exit[];
