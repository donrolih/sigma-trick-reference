##
using Plots, LaTeXStrings, DelimitedFiles
using Plots.PlotMeasures
default(lw=1.5, legendfontsize=10, fontfamily="Computer Modern")
##
G_u = readdlm("c-G-u.dat", '\t', ComplexF64)
Fl_u = readdlm("c-Fl-u.dat", '\t', ComplexF64)
I_u = readdlm("c-I-u.dat", '\t', ComplexF64)
Fr_u = readdlm("c-Fr-u.dat", '\t', ComplexF64)

awold_u = readdlm("aw-old-u.dat", '\t', ComplexF64)
sigmaold_u = readdlm("sigma-old-u.dat", '\t', ComplexF64)

aw_u = readdlm("aw-u.dat", '\t', ComplexF64)
sigma_u = readdlm("sigma-u.dat", '\t', ComplexF64)

G_d = readdlm("c-G-d.dat", '\t', ComplexF64)
Fl_d = readdlm("c-Fl-d.dat", '\t', ComplexF64)
I_d = readdlm("c-I-d.dat", '\t', ComplexF64)
Fr_d = readdlm("c-Fr-d.dat", '\t', ComplexF64)

awold_d = readdlm("aw-old-d.dat", '\t', ComplexF64)
sigmaold_d = readdlm("sigma-old-d.dat", '\t', ComplexF64)

aw_d = readdlm("aw-d.dat", '\t', ComplexF64)
sigma_d = readdlm("sigma-d.dat", '\t', ComplexF64)

ωs = real.(G_u[:, 1]);
##
# correlators
# G
p1 = plot(ωs, [(-1/pi)*real.(G_u[:, 2]) (-1/pi)*imag.(G_u[:, 2])], label=[L"\mathrm{Re}\, G_\uparrow" L"\mathrm{Im}\, G_\uparrow"])
xlabel!(p1, L"\omega")

# F
p2 = plot(ωs, [(-1/pi)*real.(Fl_u[:, 2]) (-1/pi)*imag.(Fl_u[:, 2]) (-1/pi)*real.(Fr_u[:, 2]) (-1/pi)*imag.(Fr_u[:, 2])], label=[L"\mathrm{Re}\, F^L_\uparrow" L"\mathrm{Im}\, F^L_\uparrow" L"\mathrm{Re}\, F^R_\uparrow" L"\mathrm{Im}\, F^R_\uparrow"])
xlabel!(p2, L"\omega")

p22 = plot(ωs, [(-1/pi)*real.(I_u[:, 2]) (-1/pi)*imag.(I_u[:, 2])], label=[L"\mathrm{Re}\, I_\uparrow" L"\mathrm{Im}\, I_\uparrow"])
xlabel!(p22, L"\omega")

p3 = plot(p1, p2, p22, layout=(1, 3), size=(1200, 500), margin=3mm)
##
p4 = plot(ωs, [real.(sigma_u[:, 2]) imag.(sigma_u[:, 2]) real.(sigmaold_u[:, 2]) imag.(sigmaold_u[:, 2])], label=[L"\mathrm{Re}\,\Sigma_\uparrow" L"\mathrm{Im}\,\Sigma_\uparrow" L"\mathrm{Re}\,\Sigma_\uparrow^\mathrm{old}" L"\mathrm{Im}\,\Sigma_\uparrow^\mathrm{old}"])
xlabel!(p4, L"\omega")
##
p5 = plot(ωs, [real.(aw_u[:, 2]) real.(awold_u[:, 2])], label=[L"spectral function $A_\uparrow(\omega)$"  L"spectral function $A^\mathrm{old}_\uparrow(\omega)$"], legend=:topleft)
xlabel!(p5, L"\omega")

##
# correlators
# G
p6 = plot(ωs, [(-1/pi)*real.(G_d[:, 2]) (-1/pi)*imag.(G_d[:, 2])], label=[L"\mathrm{Re}\, G_\downarrow" L"\mathrm{Im}\, G_\downarrow"])
xlabel!(p6, L"\omega")

# F
p7 = plot(ωs, [(-1/pi)*real.(Fl_d[:, 2]) (-1/pi)*imag.(Fl_d[:, 2]) (-1/pi)*real.(Fr_d[:, 2]) (-1/pi)*imag.(Fr_d[:, 2])], label=[L"\mathrm{Re}\, F^L_\downarrow" L"\mathrm{Im}\, F^L_\downarrow" L"\mathrm{Re}\, F^R_\downarrow" L"\mathrm{Im}\, F^R_\downarrow"])
xlabel!(p7, L"\omega")

p77 = plot(ωs, [(-1/pi)*real.(I_d[:, 2]) (-1/pi)*imag.(I_d[:, 2])], label=[L"\mathrm{Re}\, I_\downarrow" L"\mathrm{Im}\, I_\downarrow"])
xlabel!(p77, L"\omega")

p8 = plot(p6, p7, p77, layout=(1, 3), size=(1200, 500), margin=3mm)
##
p9 = plot(ωs, [real.(sigma_d[:, 2]) imag.(sigma_d[:, 2]) real.(sigmaold_d[:, 2]) imag.(sigmaold_d[:, 2])], label=[L"\mathrm{Re}\,\Sigma_\downarrow" L"\mathrm{Im}\,\Sigma_\downarrow" L"\mathrm{Re}\,\Sigma_\downarrow^\mathrm{old}" L"\mathrm{Im}\,\Sigma_\downarrow^\mathrm{old}"])
xlabel!(p9, L"\omega")
##
p10 = plot(ωs, [real.(aw_d[:, 2]) real.(awold_d[:, 2])], label=[L"spectral function $A_\downarrow(\omega)$"  L"spectral function $A^\mathrm{old}_\downarrow(\omega)$"], legend=:topleft)
xlabel!(p10, L"\omega")
##
plots = [p3, p4, p5, p8, p9, p10]
names = ["fig-correlators-u.pdf", "fig-sigma-u.pdf", "fig-aw-u.pdf",
         "fig-correlators-d.pdf", "fig-sigma-d.pdf", "fig-aw-d.pdf"
]

map(savefig, plots, names);