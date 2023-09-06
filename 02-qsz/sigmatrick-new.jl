# Calculate the self-energy and the spectral function via the self-energy trick
# This is the new approach by Kugler 
# QSZ symmetry type (scalar case)

using DelimitedFiles
using LinearAlgebra

###############
# Function definitions
###############

function readhybri()
    Δ11 = readdlm("Delta.dat")
    Δ11re = readdlm("Delta-re.dat")

    rows = size(Δ11, 1)
    Δ = zeros(ComplexF64, rows, 2)
    Δ[:, 1] = Δ11[:, 1]
    Δ[:, 2] = Δ11re[:, 2] .+ im*Δ11[:, 2]
    return Δ
end

###############

Nz = readchomp(`getNz`)
path = "spectral_Nz" * Nz

Δ = readhybri()

sigmaH = readdlm("custom.avg", skipstart=2)[:, 4:5]

sigmaHdict = Dict("u" => sigmaH[1, 2], "d" => sigmaH[1, 1])

for spin in ["u", "d"]
    println("Performing sigma trick for spin species $(spin)!")
    prefix = "spec_FDM_dens"
    suffix = ".dat"

    ωs = Δ[:, 1]

    imG = readdlm(joinpath(path, prefix * "_" * "A_d-A_d" * "-" * spin * suffix))
    reG = readdlm(joinpath(path, prefix * "_" * "A_d-A_d" * "-" * spin * "-re" * suffix))

    imFl = readdlm(joinpath(path, prefix * "_" * "self_d-A_d" * "-" * spin * suffix))
    reFl = readdlm(joinpath(path, prefix * "_" * "self_d-A_d" * "-" * spin * "-re" * suffix))

    imFr = readdlm(joinpath(path, prefix * "_" * "A_d-self_d" * "-" * spin * suffix))
    reFr = readdlm(joinpath(path, prefix * "_" * "A_d-self_d" * "-" * spin * "-re" * suffix))

    imI = readdlm(joinpath(path, prefix * "_" * "self_d-self_d" * "-" * spin * suffix))
    reI = readdlm(joinpath(path, prefix * "_" * "self_d-self_d" * "-" * spin * "-re" * suffix))

    G = -pi*(reG[:, 2] .+ im*imG[:, 2])
    Fl = -pi*(reFl[:, 2] .+ im*imFl[:, 2])
    Fr = -pi*(reFr[:, 2] .+ im*imFr[:, 2])
    I = -pi*(reI[:, 2] .+ im*imI[:, 2])

    open("c-G-" * spin * ".dat", "w") do io
        corr = [complex(ωs) G]
        writedlm(io, corr)
    end

    open("c-Fl-" * spin * ".dat", "w") do io
        corr = [complex(ωs) Fl]
        writedlm(io, corr)
    end

    open("c-Fr-" * spin * ".dat", "w") do io
        corr = [complex(ωs) Fr]
        writedlm(io, corr)
    end

    open("c-I-" * spin * ".dat", "w") do io
        corr = [complex(ωs) I]
        writedlm(io, corr)
    end
    

    # Calculate the self-energy

    sigmaold = Fl ./ G 
    sigma = sigmaHdict[spin] .+ I .- (Fl .^ 2) ./ G 

    # Green's function
    gfold = 1. ./ (ωs .+ Δ[:, 2] .- sigmaold)
    # Imaginary part of the Green's function (spectral function)
    awold = -1/pi * imag.(gfold)

    # Green's function
    gf = 1. ./ (ωs .+ Δ[:, 2] .- sigma)
    # Imaginary part of the Green's function (spectral function)
    aw = -1/pi * imag.(gf)


    open("sigma-old-" * spin * ".dat", "w") do io
        writedlm(io, [ωs sigmaold])
    end

    open("aw-old-" * spin * ".dat", "w") do io
        writedlm(io, [ωs awold])
    end

    open("sigma-" * spin * ".dat", "w") do io
        writedlm(io, [ωs sigma])
    end

    open("aw-" * spin * ".dat", "w") do io
        writedlm(io, [ωs aw])
    end

end