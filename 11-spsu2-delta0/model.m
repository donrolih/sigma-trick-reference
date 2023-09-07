def1ch[1];

H = H0 + Hc + H1;

(* All operators which contain d[], except hybridization (Hc). *)
Hselfd = H1;

(* These operators are necessary for second order Green's functions *)
selfopd =  ( Chop @ Expand @ komutator[Hselfd /. params, d[#1, #2]] )&;
selfopcd = ( Chop @ Expand @ komutator[Hselfd /. params, ((-1)^#2 d[1-#1, 1-#2])] )&;

(* selfopd =  ( Chop @ Expand @ komutator[ d[#1, #2], Hselfd /. params] )&;
selfopcd = ( Chop @ Expand @ komutator[((-1)^#2 d[1-#1, 1-#2]), Hselfd /. params] )&; *)

(* This is the constant Hartree term needed for the self-energy *)
SigmaHd = Expand @ antikomutator[ selfopd[#1, #2], d[#3, #4] ] /. params &;
SigmaHdAvg11 := Expand @ (SigmaHd[CR, UP, AN, UP]+SigmaHd[CR, DO, AN, DO])/2;
SigmaHdAvg12 := Expand @ (SigmaHd[AN, DO, AN, UP]+SigmaHd[AN, UP, AN, DO])/2;
SigmaHdAvg21 := Expand @ (SigmaHd[CR, UP, CR, DO]+SigmaHd[CR, DO, CR, UP])/2;
SigmaHdAvg22 := Expand @ (SigmaHd[AN, DO, CR, DO]+SigmaHd[AN, UP, CR, UP])/2;


Print["SigmaH[AN, UP, CR, UP]=", SigmaHd[CR, UP, AN, UP] ];
Print["SigmaH[AN, UP, AN, DO]=", SigmaHd[AN, DO, AN, UP] ];
Print["SigmaH[CR, DO, CR, UP]=", SigmaHd[CR, UP, CR, DO] ];
Print["SigmaH[CR, DO, AN, DO]=", SigmaHd[AN, DO, CR, DO] ];

Print["SigmaHdAvg11=", SigmaHdAvg11];
Print["SigmaHdAvg12=", SigmaHdAvg12];
Print["SigmaHdAvg21=", SigmaHdAvg21];
Print["SigmaHdAvg22=", SigmaHdAvg22];

(* Evaluate *)
Print["selfopd[CR,UP]=", selfopd[CR, UP]];
Print["selfopd[CR,DO]=", selfopd[CR, DO]];
Print["selfopd[AN,UP]=", selfopd[AN, UP]];
Print["selfopd[AN,DO]=", selfopd[AN, DO]];

Print["selfopcd[CR,UP]=", selfopcd[CR, UP]];
Print["selfopcd[CR,DO]=", selfopcd[CR, DO]];
Print["selfopcd[AN,UP]=", selfopcd[AN, UP]];
Print["selfopcd[AN,DO]=", selfopcd[AN, DO]];
