def1ch[1];

H1 += B spinz[d[]];

H = H0 + Hc + H1;

(* All operators which contain d[], except hybridization (Hc). *)
Hselfd = H1;

selfopd =  ( Chop @ Expand @ komutator[Hselfd /. params, d[#1, #2]] )&;
selfopcd = ( Chop @ Expand @ komutator[Hselfd /. params, ((-1)^#2 d[1-#1, 1-#2])] )&;

SigmaHd = Expand @ antikomutator[ selfopd[#1, #2], d[#3, #4] ] /. params &;

Print["SigmaH[CR, UP, AN, UP]=", SigmaHd[CR, UP, AN, UP] ];
Print["SigmaH[AN, DO, AN, UP]=", SigmaHd[AN, DO, AN, UP] ];
Print["SigmaH[CR, UP, CR, DO]=", SigmaHd[CR, UP, CR, DO] ];
Print["SigmaH[AN, DO, CR, DO]=", SigmaHd[AN, DO, CR, DO] ];

Print["SigmaH[CR, DO, AN, DO]=", SigmaHd[CR, DO, AN, DO] ];
Print["SigmaH[AN, UP, AN, DO]=", SigmaHd[AN, UP, AN, DO] ];
Print["SigmaH[CR, DO, CR, UP]=", SigmaHd[CR, DO, CR, UP] ];
Print["SigmaH[AN, UP, CR, UP]=", SigmaHd[AN, UP, CR, UP] ];

(* Evaluate *)
Print["selfopd[CR,UP]=", selfopd[CR, UP]];
Print["selfopd[CR,DO]=", selfopd[CR, DO]];
Print["selfopd[AN,UP]=", selfopd[AN, UP]];
Print["selfopd[AN,DO]=", selfopd[AN, DO]];

Print["selfopcd[CR,UP]=", selfopcd[CR, UP]];
Print["selfopcd[CR,DO]=", selfopcd[CR, DO]];
Print["selfopcd[AN,UP]=", selfopcd[AN, UP]];
Print["selfopcd[AN,DO]=", selfopcd[AN, DO]];
