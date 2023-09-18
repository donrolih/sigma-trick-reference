using PythonPlot, LaTeXStrings
using DelimitedFiles
pyplot.style.use(["don-custom"]);
##
custom = readdlm("custom.avg", skipstart=2)
sigmaH = custom[3:6]
# it should be like this
sigmaHm = [sigmaH[1] sigmaH[3]; sigmaH[2] sigmaH[4]];
sigmaHC = custom[7:10]
# it should be like this
sigmaHCm = [sigmaHC[1] sigmaHC[3]; sigmaHC[2] sigmaHC[4]];
##
G_u = readdlm("c-G-u.dat", '\t', ComplexF64; skipstart=1)
Fl_u = readdlm("c-Fl-u.dat", '\t', ComplexF64; skipstart=1)
I_u = readdlm("c-I-u.dat", '\t', ComplexF64; skipstart=1)
Fr_u = readdlm("c-Fr-u.dat", '\t', ComplexF64; skipstart=1)

awold_u = readdlm("aw-old-u.dat", Float64)
sigmaold_u = readdlm("sigma-old-u.dat", '\t', ComplexF64)

aw_u = readdlm("aw-u.dat", Float64)
sigma_u = readdlm("sigma-u.dat", '\t', ComplexF64)

G_d = readdlm("c-G-d.dat", '\t', ComplexF64; skipstart=1)
Fl_d = readdlm("c-Fl-d.dat", '\t', ComplexF64; skipstart=1)
I_d = readdlm("c-I-d.dat", '\t', ComplexF64; skipstart=1)
Fr_d = readdlm("c-Fr-d.dat", '\t', ComplexF64; skipstart=1)

awold_d = readdlm("aw-old-d.dat", Float64)
sigmaold_d = readdlm("sigma-old-d.dat", '\t', ComplexF64)

aw_d = readdlm("aw-d.dat", Float64)
sigma_d = readdlm("sigma-d.dat", '\t', ComplexF64)

ωs = real.(G_u[:, 1]);

##
fig, ax = pyplot.subplots()

ax.plot(ωs, G_u[:, 5], label="u")
ax.plot(ωs, G_d[:, 5], label="d")

display(fig)
##
# Spin UP
fig, ax = pyplot.subplots(ncols=2, figsize=(9,5))

label_G11 = L"A_{11}= -\frac{1}{\pi}\mathrm{Im} \langle\langle d_\uparrow; d_\uparrow^\dagger \rangle\rangle_z"

label_G31 = L"A_{31}= -\frac{1}{\pi}\mathrm{Im}\langle\langle -d_\downarrow^\dagger; d_\uparrow^\dagger \rangle\rangle_z"

label_G13 = L"A_{13}= -\frac{1}{\pi}\mathrm{Im}\langle\langle d_\uparrow; -d_\downarrow \rangle\rangle_z"

label_G33 = L"A_{33}= -\frac{1}{\pi}\mathrm{Im}\langle\langle d_\downarrow^\dagger; d_\downarrow \rangle\rangle_z"

ax[0].plot(ωs, aw_u[:, 2], label=label_G11)
ax[0].plot(ωs, aw_u[:, 5], label=label_G33)
ax[1].plot(ωs, aw_u[:, 3], label=label_G13)
ax[1].plot(ωs, aw_u[:, 4], label=label_G31)

# ax.plot(ωs, aw_u[:, 2], label="new up")
ax[0].set_xlim(-0.3, 0.3)
ax[1].set_xlim(-0.003, 0.003)
# ax.plot(ωs, awold[:, 5])
# ax.set_xlim(-0.00025, 0.00025)
ax[0].legend(loc="upper left")
ax[1].legend()
fig.tight_layout()
display(fig)
savefig("fig-aw-u.pdf")
##
# Spin DOWN
fig, ax = pyplot.subplots(ncols=2, figsize=(9,5))

label_G22 = L"A_{22}= -\frac{1}{\pi}\mathrm{Im}\langle\langle d_\downarrow; d_\downarrow^\dagger \rangle\rangle_z"

label_G42 = L"A_{42}= -\frac{1}{\pi}\mathrm{Im}\langle\langle d_\uparrow^\dagger; d_\downarrow^\dagger \rangle\rangle_z"

label_G24 = L"A_{24}= -\frac{1}{\pi}\mathrm{Im}\langle\langle d_\downarrow; d_\uparrow \rangle\rangle_z"

label_G44 = L"A_{44}= -\frac{1}{\pi}\mathrm{Im}\langle\langle d_\uparrow^\dagger; d_\uparrow \rangle\rangle_z"

ax[0].plot(ωs, aw_d[:, 2], label=label_G22)
ax[0].plot(ωs, aw_d[:, 5], label=label_G44)
ax[1].plot(ωs, aw_d[:, 3], label=label_G24)
ax[1].plot(ωs, aw_d[:, 4], label=label_G42)

# ax.plot(ωs, aw_u[:, 2], label="new up")
ax[0].set_xlim(-0.3, 0.3)
ax[1].set_xlim(-0.003, 0.003)
# ax.plot(ωs, awold[:, 5])
# ax.set_xlim(-0.00025, 0.00025)
ax[0].legend(loc="upper left")
ax[1].legend()
fig.tight_layout()
display(fig)
savefig("fig-aw-d.pdf")
##