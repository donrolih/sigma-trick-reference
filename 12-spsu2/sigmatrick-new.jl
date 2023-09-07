# Calculate the self-energy and the spectral function via the self-energy trick
# This is the new approach by Kugler 
# SPSU2 symmetry type (2x2 matrix case)

# The correlators, self-energies and spectral functions
# are stored into .dat files with 5 columns
# the first column are the frequencies
# and then the elements of the matrices in the
# order: 11, 21, 12, 22.
# This is because Julia is a column major language
# and therefore reshape([a11, a21, a12, a22], 2, 2) returns
# a correct matrix [a11 a12; a21 a22]

using DelimitedFiles
using LinearAlgebra

###############
# Function definitions
###############

function readhybri()
    Δ11 = readdlm("Delta11.dat")
    Δ12 = readdlm("Delta12.dat")
    Δ21 = readdlm("Delta21.dat")
    Δ22 = readdlm("Delta22.dat")

    Δ11re = readdlm("Delta11-re.dat")
    Δ12re = readdlm("Delta12-re.dat")
    Δ21re = readdlm("Delta21-re.dat")
    Δ22re = readdlm("Delta22-re.dat")

    rows = size(Δ11, 1)
    Δ = zeros(ComplexF64, rows, 5)
    Δ[:, 1] = Δ11[:, 1]
    Δ[:, 2] = Δ11re[:, 2] .+ im*Δ11[:, 2]
    Δ[:, 4] = Δ12re[:, 2] .+ im*Δ12[:, 2]
    Δ[:, 3] = Δ21re[:, 2] .+ im*Δ21[:, 2]
    Δ[:, 5] = Δ22re[:, 2] .+ im*Δ22[:, 2]

    return Δ
end

function makespecname(opsc::Vector{String}, ops::Vector{String})
    n = size(opsc, 1)
    m = size(ops, 1)

    specmatrix = fill("", n, m)

    for (j, op) in enumerate(ops)
        for (i, opc) in enumerate(opsc)
            specmatrix[i, j] = opc * "-" * op
        end
    end

    return specmatrix
end

function readspecs(imp_ops, self_ops, path; prefix="spec_FDM_dens", suffix=".dat")
    names = Dict()
    names["G"] = makespecname(imp_ops, imp_ops)
    names["Fl"] = makespecname(self_ops, imp_ops)
    names["Fr"] = makespecname(imp_ops, self_ops)
    names["I"] = makespecname(self_ops, self_ops)
    
    n = length(imp_ops) * length(imp_ops)

    specdata = Dict()

    for matrix in ["G", "Fl", "Fr", "I"]
        corr = ComplexF64[]
        for i in 1:n
            name = names[matrix][i]
            filenameim = prefix * "_" * name * suffix
            dataim = readdlm(joinpath(path, filenameim))
            filenamere = prefix * "_" * name * "-re" * suffix
            datare = readdlm(joinpath(path, filenamere))
            col = -pi*(datare[:, 2] .+ im*dataim[:, 2])
            if i == 1
                corr = [dataim[:, 1] col]
            else
                corr = hcat(corr, col)
            end
        end
        specdata[matrix] = corr
    end

    G = specdata["G"]
    Fl = specdata["Fl"]
    Fr = specdata["Fr"]
    I = specdata["I"]

    return G, Fl, Fr, I
end

function sigmatrick(Δ, sigmaHm, G, Fl, Fr, I)
    ωs = Δ[:, 1]

    sigma = zeros(ComplexF64, 1, 5)
    sigmaold = zeros(ComplexF64, 1, 5)

    aw = zeros(Float64, 1, 5)
    awold = zeros(Float64, 1, 5)

    τ₃ = [(1. + 0. * im) (0. + 0. * im); (0. + 0. * im) (-1. + 0. * im)]

    for (i, ω) in enumerate(ωs)
        G_m = reshape(G[i, 2:end], 2, 2)
        Fl_m = reshape(Fl[i, 2:end], 2, 2)
        Fr_m = reshape(Fr[i, 2:end], 2, 2)
        I_m = reshape(I[i, 2:end], 2, 2)

        detG = G_m[1, 1]*G_m[2, 2] - G_m[1, 2]*G_m[2, 1]
        G_m_inv = (1/detG) .* [G_m[2, 2] -G_m[1, 2]; -G_m[2, 1] G_m[1, 1]]

        # Kugler's method: Σ = SigmaHd .+ I .- Fl ./ G .* Fr
        Σ_m = sigmaHm + I_m - Fl_m * G_m_inv * Fr_m
        # Old sigma trick (Bulla, Hewson, Pruschke): Σ = Fl/G = G\Fr
        # Both Fl and Fr should give the same result (up to numerical error)
        Σ_old_m = Fl_m * G_m_inv

        Δ_m = reshape(Δ[i, 2:end], 2, 2)

        Gfinv = ω*LinearAlgebra.I + Δ_m - Σ_m
        detGf = Gfinv[1, 1]*Gfinv[2, 2] - Gfinv[1, 2]*Gfinv[2, 1] 
        Gf_m = (1/detGf) * [Gfinv[2, 2] -Gfinv[1, 2]; -Gfinv[2, 1] Gfinv[1, 1]]

        Gfinvold = ω*LinearAlgebra.I + Δ_m - Σ_old_m
        detGfold = Gfinvold[1, 1]*Gfinvold[2, 2] - Gfinvold[1, 2]*Gfinvold[2, 1] 
        Gf_old_m = (1/detGfold) * [Gfinvold[2, 2] -Gfinvold[1, 2]; -Gfinvold[2, 1] Gfinvold[1, 1]]

        A_m = -(1/pi)*imag.(Gf_m)
        A_old_m = -(1/pi)*imag.(Gf_old_m)
        if i == 1
            sigma = reshape(pushfirst!(vec(reshape(Σ_m, 1, 4)), ω), 1, 5)
            sigmaold = reshape(pushfirst!(vec(reshape(Σ_old_m, 1, 4)), ω), 1, 5)
            aw = reshape(pushfirst!(vec(reshape(A_m, 1, 4)), ω), 1, 5)
            awold = reshape(pushfirst!(vec(reshape(A_old_m, 1, 4)), ω), 1, 5)
        else
            sigma = vcat(sigma, reshape(pushfirst!(vec(reshape(Σ_m, 1, 4)), ω), 1, 5))
            sigmaold = vcat(sigmaold, reshape(pushfirst!(vec(reshape(Σ_old_m, 1, 4)), ω), 1, 5))
            aw = vcat(aw, reshape(pushfirst!(vec(reshape(A_m, 1, 4)), ω), 1, 5))
            awold = vcat(awold, reshape(pushfirst!(vec(reshape(A_old_m, 1, 4)), ω), 1, 5))
        end
    end
    return sigma, sigmaold, aw, awold
end

###############

imp_ops = ["A_d", "Ac_d"]
self_ops = ["self_d", "selfc_d"]

Nz = readchomp(`getNz`)
path = "spectral_Nz" * Nz

Δ = readhybri()
G, Fl, Fr, I = readspecs(imp_ops, self_ops, path)

open("c-G.dat", "w") do io
    header = reshape(makespecname(imp_ops, imp_ops), 1, 4)
    header = hcat("w", header)
    writedlm(io, [header; G])
end

open("c-Fl.dat", "w") do io
    header = reshape(makespecname(self_ops, imp_ops), 1, 4)
    header = hcat("w", header)
    writedlm(io, [header; Fl])
end

open("c-Fr.dat", "w") do io
    header = reshape(makespecname(imp_ops, self_ops), 1, 4)
    header = hcat("w", header)
    writedlm(io, [header; Fr])
end

open("c-I.dat", "w") do io
    header = reshape(makespecname(self_ops, self_ops), 1, 4)
    header = hcat("w", header)
    writedlm(io, [header; I])
end

custom = readdlm("custom.avg", skipstart=2)

# careful!
# element 11 corresponds to up - up spin species
# element 12 to up - down
# element 21 to down - up
# element 22 to down - down
# matrix form
sigmaH11 = parse(Float64, readchomp(`extractcolumn custom.avg SigmaHdAvg11`))
sigmaH12 = parse(Float64, readchomp(`extractcolumn custom.avg SigmaHdAvg12`))
sigmaH21 = parse(Float64, readchomp(`extractcolumn custom.avg SigmaHdAvg21`))
sigmaH22 = parse(Float64, readchomp(`extractcolumn custom.avg SigmaHdAvg22`))

##
# it should be like this
sigmaHm = [sigmaH11 sigmaH12; sigmaH21 sigmaH22];

sigma, sigmaold, aw, awold = sigmatrick(Δ, sigmaHm, G, Fl, Fr, I)

open("sigma.dat", "w") do io
    writedlm(io, sigma)
end

open("sigma-old.dat", "w") do io
    writedlm(io, sigmaold)
end

open("aw.dat", "w") do io
    writedlm(io, aw)
end

open("aw-old.dat", "w") do io
    writedlm(io, awold)
end