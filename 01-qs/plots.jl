##
using Plots, LaTeXStrings, DelimitedFiles
using Plots.PlotMeasures
default(lw=1.5, legendfontsize=10, fontfamily="Computer Modern")
##
G = readdlm("c-G.dat", '\t', ComplexF64)
Fl = readdlm("c-Fl.dat", '\t', ComplexF64)
I = readdlm("c-I.dat", '\t', ComplexF64)
Fr = readdlm("c-Fr.dat", '\t', ComplexF64)

awold = readdlm("aw-old.dat", '\t', ComplexF64)
sigmaold = readdlm("sigma-old.dat", '\t', ComplexF64)

aw = readdlm("aw.dat", '\t', ComplexF64)
sigma = readdlm("sigma.dat", '\t', ComplexF64)

ωs = real.(G[:, 1]);
##
# correlators
# G
p1 = plot(ωs, [(-1/pi)*real.(G[:, 2]) (-1/pi)*imag.(G[:, 2])], label=[L"\mathrm{Re}\, G" L"\mathrm{Im}\, G"])
xlabel!(p1, L"\omega")

# F
p2 = plot(ωs, [(-1/pi)*real.(Fl[:, 2]) (-1/pi)*imag.(Fl[:, 2]) (-1/pi)*real.(Fr[:, 2]) (-1/pi)*imag.(Fr[:, 2])], label=[L"\mathrm{Re}\, F^L" L"\mathrm{Im}\, F^L" L"\mathrm{Re}\, F^R" L"\mathrm{Im}\, F^R"])
xlabel!(p2, L"\omega")

p22 = plot(ωs, [(-1/pi)*real.(I[:, 2]) (-1/pi)*imag.(I[:, 2])], label=[L"\mathrm{Re}\, I" L"\mathrm{Im}\, I"])
xlabel!(p22, L"\omega")

p3 = plot(p1, p2, p22, layout=(1, 3), size=(1200, 500), margin=3mm)
##
p4 = plot(ωs, [real.(sigma[:, 2]) imag.(sigma[:, 2]) real.(sigmaold[:, 2]) imag.(sigmaold[:, 2])], label=[L"\mathrm{Re}\,\Sigma" L"\mathrm{Im}\,\Sigma" L"\mathrm{Re}\,\Sigma_\mathrm{old}" L"\mathrm{Im}\,\Sigma_\mathrm{old}"])
xlabel!(p4, L"\omega")
##
p5 = plot(ωs, [real.(aw[:, 2]) real.(awold[:, 2])], label=[L"spectral function $A(\omega)$"  L"spectral function $A_\mathrm{old}(\omega)$"])
xlabel!(p5, L"\omega")

##
plots = [p3, p4, p5]
names = ["fig-correlators.pdf", "fig-sigma.pdf", "fig-aw.pdf"]

map(savefig, plots, names);