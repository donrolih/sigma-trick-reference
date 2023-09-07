#!/usr/bin/env perl
use strict;
use warnings;

my $Gamma = `getparam Gamma param.loop`;
my $bcsgap = `getparam bcsgap param.loop`;
chomp($Gamma, $bcsgap);

my $scr = "script.tmp";

open (F, ">$scr") or die;
print F <<EOF;

Delta = $bcsgap;
gamma = $Gamma;

(* finite flat band *)

(* NOTE THE BCS GAP PHASE CONVENTION, i.e., THE MINUS SIGN! *)
h[z_] := {{-((z ArcCot[Sqrt[-z^2 + Delta^2]])/Sqrt[-z^2 + Delta^2]), 
(-Delta ArcCot[Sqrt[-z^2 + Delta^2]])/Sqrt[-z^2 + Delta^2]}, 
{(-Delta ArcCot[Sqrt[-z^2 + Delta^2]])/Sqrt[-z^2 + Delta^2], 
-((z ArcCot[Sqrt[-z^2 + Delta^2]])/Sqrt[-z^2 + Delta^2])}};

l = Import["Delta11.dat.flat", "Table"];
f = Map[{First[#], -(2/Pi) gamma Im[ h[First[#]+10^-10 I][[1,1]] ]}&, l];
Export["Delta11.dat", f];

l = Import["Delta12.dat.flat", "Table"];
f = Map[{First[#], -(2/Pi) gamma Im[ h[First[#]+10^-10 I][[1,2]] ]}&, l];
Export["Delta12.dat", f];

EOF

system "math -batchinput -batchoutput <$scr >/dev/null 2>/dev/null";

