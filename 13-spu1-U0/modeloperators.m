Module[{t},
  t = {};
  t = Join[t, mtDoubletOp["Ac_d", ((-1)^#2 d[1-#1, 1-#2])&] ];
  MPVCFAST = False;
  t = Join[t, mtDoubletOp["selfc_d", selfopcd ]];
  t = Join[t, mtSingletOp["SigmaHd11-u", SigmaHd[CR, UP, AN, UP] ]];
  t = Join[t, mtSingletOp["SigmaHd12-u", SigmaHd[AN, DO, AN, UP] ]];
  t = Join[t, mtSingletOp["SigmaHd21-u", SigmaHd[CR, UP, CR, DO] ]];
  t = Join[t, mtSingletOp["SigmaHd22-u", SigmaHd[AN, DO, CR, DO] ]];
  t = Join[t, mtSingletOp["SigmaHd11-d", SigmaHd[CR, DO, AN, DO] ]];
  t = Join[t, mtSingletOp["SigmaHd12-d", SigmaHd[AN, UP, AN, DO] ]];
  t = Join[t, mtSingletOp["SigmaHd21-d", SigmaHd[CR, DO, CR, UP] ]];
  t = Join[t, mtSingletOp["SigmaHd22-d", SigmaHd[AN, UP, CR, UP] ]];
  MPVCFAST = True;
  texportable = t;
];
texportable
