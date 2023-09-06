Module[{t},
  t = {};
  MPVCFAST = False;
  t = Join[t, mtSingletOp["SigmaHd",   SigmaHdAvg ]];
  t = Join[t, mtSingletOp["SigmaHd-u-u", SigmaHd[UP, UP] ]];
  t = Join[t, mtSingletOp["SigmaHd-u-d", SigmaHd[UP, DO] ]];
  t = Join[t, mtSingletOp["SigmaHd-d-u", SigmaHd[DO, UP] ]];
  t = Join[t, mtSingletOp["SigmaHd-d-d", SigmaHd[DO, DO] ]];
  MPVCFAST = True;
  texportable = t;
];
texportable
