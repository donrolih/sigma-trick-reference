##
using Plots, LaTeXStrings, DelimitedFiles
using Plots.PlotMeasures
default(lw=1.5, legendfontsize=10, fontfamily="Computer Modern")
##
G = readdlm("c-G.dat", '\t', ComplexF64; skipstart=1)
Fl = readdlm("c-Fl.dat", '\t', ComplexF64; skipstart=1)
I = readdlm("c-I.dat", '\t', ComplexF64; skipstart=1)
Fr = readdlm("c-Fr.dat", '\t', ComplexF64; skipstart=1)

awold = readdlm("aw-old.dat", Float64)
sigmaold = readdlm("sigma-old.dat", '\t', ComplexF64)

aw = readdlm("aw.dat", Float64)
sigma = readdlm("sigma.dat", '\t', ComplexF64)

ωs = real.(G[:, 1]);
##
# correlators
# G
p1 = plot(ωs, (-1/pi)*real.(G[:, 2:end]), label=[L"\mathrm{Re}\, G_{11}" L"\mathrm{Re}\, G_{21}" L"\mathrm{Re}\, G_{12}" L"\mathrm{Re}\, G_{22}"])
xlabel!(p1, L"\omega")

p2 = plot(ωs, (-1/pi)*imag.(G[:, 2:end]), label=[L"\mathrm{Im}\, G_{11}" L"\mathrm{Im}\, G_{21}" L"\mathrm{Im}\, G_{12}" L"\mathrm{Im}\, G_{22}"])
xlabel!(p2, L"\omega")

p3 = plot(p1, p2, layout=(1, 2), size=(1000, 400), margin=3mm)
##
# F
p4 = plot(ωs, (-1/pi)*real.(Fl[:, 2:end]), label=[L"\mathrm{Re}\, F^{\mathrm{L}}_{11}" L"\mathrm{Re}\, F^{\mathrm{L}}_{21}" L"\mathrm{Re}\, F^{\mathrm{L}}_{12}" L"\mathrm{Re}\, F^{\mathrm{L}}_{22}"])
xlabel!(p4, L"\omega")

p5 = plot(ωs, (-1/pi)*imag.(Fl[:, 2:end]), label=[L"\mathrm{Im}\, F^{\mathrm{L}}_{11}" L"\mathrm{Im}\, F^{\mathrm{L}}_{21}" L"\mathrm{Im}\, F^{\mathrm{L}}_{12}" L"\mathrm{Im}\, F^{\mathrm{L}}_{22}"])
xlabel!(p5, L"\omega")

p6 = plot(p4, p5, layout=(1, 2), size=(1000, 400), margin=3mm)
##
#I
p7 = plot(ωs, (-1/pi)*real.(I[:, 2:end]), label=[L"\mathrm{Re}\, I_{11}" L"\mathrm{Re}\, I_{21}" L"\mathrm{Re}\, I_{12}" L"\mathrm{Re}\, I_{22}"])
xlabel!(p7, L"\omega")

p8 = plot(ωs, (-1/pi)*imag.(I[:, 2:end]), label=[L"\mathrm{Im}\,I_{11}" L"\mathrm{Im}\, I_{21}" L"\mathrm{Im}\, I_{12}" L"\mathrm{Im}\, I_{22}"])
xlabel!(p8, L"\omega")

p9 = plot(p7, p8, layout=(1, 2), size=(1000, 400), margin=3mm)
##

labels=[L"\mathrm{Re}\, \Sigma_{11}" L"\mathrm{Re}\, \Sigma_{21}" L"\mathrm{Re}\, \Sigma_{12}" L"\mathrm{Re}\, \Sigma_{22}"]
p10 = plot(ωs, real.(sigma[:, 2:end]), label=labels)
xlabel!(p10, L"\omega")

labelsold=[L"\mathrm{Re}\, \Sigma^\mathrm{old}_{11}" L"\mathrm{Re}\, \Sigma^\mathrm{old}_{21}" L"\mathrm{Re}\, \Sigma^\mathrm{old}_{12}" L"\mathrm{Re}\, \Sigma^\mathrm{old}_{22}"]
p11 = plot(ωs, real.(sigmaold[:, 2:end]), label=labelsold)
xlabel!(p11, L"\omega")

p12 = plot(p10, p11, layout=(1, 2), size=(1000, 400), margin=3mm)
##
labels=[L"\mathrm{Im}\, \Sigma_{11}" L"\mathrm{Im}\, \Sigma_{21}" L"\mathrm{Im}\, \Sigma_{12}" L"\mathrm{Im}\, \Sigma_{22}"]
p13 = plot(ωs, imag.(sigma[:, 2:end]), label=labels)
xlabel!(p13, L"\omega")

labelsold=[L"\mathrm{Im}\, \Sigma^\mathrm{old}_{11}" L"\mathrm{Im}\, \Sigma^\mathrm{old}_{21}" L"\mathrm{Im}\, \Sigma^\mathrm{old}_{12}" L"\mathrm{Im}\, \Sigma^\mathrm{old}_{22}"]
p14 = plot(ωs, imag.(sigmaold[:, 2:end]), label=labelsold)
xlabel!(p14, L"\omega")

p15 = plot(p13, p14, layout=(1, 2), size=(1000, 400), margin=3mm)

##
labels = [L"A_{\uparrow\uparrow}" L"A_{\downarrow\uparrow}" L"A_{\uparrow\downarrow}" L"A_{\downarrow\downarrow}"]

p16 = plot(ωs, aw[:, 2:end], label=labels, xlabel=L"\omega")

labelsold = [L"A^\mathrm{old}_{\uparrow\uparrow}" L"A^\mathrm{old}_{\downarrow\uparrow}" L"A^\mathrm{old}_{\uparrow\downarrow}" L"A^\mathrm{old}_{\downarrow\downarrow}"]

p17 = plot(ωs, awold[:, 2:end], label=labelsold, xlabel=L"\omega")

p18 = plot(p16, p17, layout=(1, 2), size=(1000, 400))
##
refsigma11re = readdlm("referencesigma/resigma11.dat", Float32)

plot(refsigma11re[:, 1], refsigma11re[:, 2])
plot!(ωs, real.(sigma[:, 2]))
plot!(ωs, real.(sigmaold[:, 2]))
##
plots = [p3, p6, p9, p12, p15, p18]
names = ["fig-correlator-G.pdf",
         "fig-correlator-F.pdf",
         "fig-correlator-I.pdf",
         "fig-sigmare.pdf",
         "fig-sigmaim.pdf",
         "fig-aw.pdf"]

map(savefig, plots, names);