using PythonPlot, LaTeXStrings
using DelimitedFiles
pyplot.style.use(["don-custom"]);
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
fig, ax = pyplot.subplots()

ax.plot(ωs, aw[:, 2])
ax.plot(ωs, awold[:, 2])
# ax.plot(ωs, awold[:, 5])
# ax.set_ylim(-0.02, 5)
display(fig)
##
fig, ax = pyplot.subplots()

ax.plot(ωs, aw[:, 5])
ax.plot(ωs, awold[:, 5])
# ax.plot(ωs, awold[:, 5])

display(fig)
##
fig, ax = pyplot.subplots(ncols=2, sharey=true, figsize=(9,4))

ax[0].plot(ωs, aw[:, 2], label=L"A_{11}(\omega)")
ax[0].plot(ωs, aw[:, 4], label=L"A_{21}(\omega)")
ax[0].plot(ωs, aw[:, 5], label=L"A_{22}(\omega)")

ax[1].plot(ωs, awold[:, 2], label=L"A_{11}^\mathrm{old}(\omega)")
ax[1].plot(ωs, awold[:, 4], label=L"A_{21}^\mathrm{old}(\omega)")
ax[1].plot(ωs, awold[:, 5], label=L"A_{22}^\mathrm{old}(\omega)")

ax[0].set_xlim(-0.6, 0.6)
ax[1].set_xlim(-0.6, 0.6)

ax[0].legend()
ax[1].legend()

fig.tight_layout(w_pad = 2)

fig.suptitle(y=1.02, L"U=0.5,\; \delta = 0.1, \; \Delta = 0.02")

display(fig)

# fig.savefig("fig-aw.pdf");
##
fig, ax = pyplot.subplots(ncols=2, sharey=true, figsize=(9,4))

ax[0].plot(ωs, sigma[:, 2], label=L"\mathrm{Re}\,\Sigma_{11}(\omega)")
ax[0].plot(ωs, sigma[:, 4], label=L"\mathrm{Re}\,\Sigma_{21}(\omega)")
ax[0].plot(ωs, sigma[:, 5], label=L"\mathrm{Re}\,\Sigma_{22}(\omega)")

ax[1].plot(ωs, sigmaold[:, 2], label=L"\mathrm{Re}\,\Sigma_{11}^\mathrm{old}(\omega)")
ax[1].plot(ωs, sigmaold[:, 4], label=L"\mathrm{Re}\,\Sigma_{21}^\mathrm{old}(\omega)")
ax[1].plot(ωs, sigmaold[:, 5], label=L"\mathrm{Re}\,\Sigma_{22}^\mathrm{old}(\omega)")

ax[0].set_xlim(-0.6, 0.6)
ax[1].set_xlim(-0.6, 0.6)

ax[0].legend()
ax[1].legend()

fig.tight_layout(w_pad = 2)

fig.suptitle(y=1.02, L"U=0.5,\; \delta = 0.1, \; \Delta = 0.02")

display(fig)

fig.savefig("fig-sigma.pdf");