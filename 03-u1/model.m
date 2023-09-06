def1ch[1];

H1 += B spinz[d[]] + Bx spinx[d[]];

H = H0 + Hc + H1;

(* All operators which contain d[], except hybridization (Hc). *)
Hselfd = H1;

selfopd = ( Chop @ Expand @ komutator[Hselfd /. params, d[#1, #2]] )&;

Print["self[CR,UP]=", selfopd[CR,UP] ];

SigmaHd = Expand @ antikomutator[ selfopd[CR, #1], d[AN, #2] ] /. params &;
SigmaHdAvg := Expand @ (SigmaHd[UP]+SigmaHd[DO])/2;

Print["SigmaH[UP, UP]=", SigmaHd[UP, UP] ];
Print["SigmaH[UP, DO]=", SigmaHd[UP, DO] ];
Print["SigmaH[DO, UP]=", SigmaHd[DO, UP] ];
Print["SigmaH[DO, DO]=", SigmaHd[DO, DO] ];
Print["SigmaH=", SigmaHdAvg ];

(* Evaluate *)
Print["selfopd[CR,UP]=", selfopd[CR, UP]];
Print["selfopd[CR,DO]=", selfopd[CR, DO]];
Print["selfopd[AN,UP]=", selfopd[AN, UP]];
Print["selfopd[AN,DO]=", selfopd[AN, DO]];
