Module[{t},
  t = {};
  t = Join[t, mtDoubletOp["Ac_d", ((-1)^#2 d[1-#1, 1-#2])&] ];
  MPVCFAST = False;
  t = Join[t, mtDoubletOp["selfc_d", selfopcd ]];
  t = Join[t, mtSingletOp["SigmaHdAvg11",   SigmaHdAvg11 ]];
  t = Join[t, mtSingletOp["SigmaHdAvg22",   SigmaHdAvg22 ]];
  t = Join[t, mtSingletOp["SigmaHdAvg12",   SigmaHdAvg12 ]];
  t = Join[t, mtSingletOp["SigmaHdAvg21",   SigmaHdAvg21 ]];
  t = Join[t, mtSingletOp["SigmaHd11", SigmaHd[CR, UP, AN, UP] ]];
  t = Join[t, mtSingletOp["SigmaHd12", SigmaHd[AN, DO, AN, UP] ]];
  t = Join[t, mtSingletOp["SigmaHd21", SigmaHd[CR, UP, CR, DO] ]];
  t = Join[t, mtSingletOp["SigmaHd22", SigmaHd[AN, DO, CR, DO] ]];
  MPVCFAST = True;
  texportable = t;
];
texportable
