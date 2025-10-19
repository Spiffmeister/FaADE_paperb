#=
    Run the MMS for the multiblock NIMROD case.
=#
using LinearAlgebra
using DataFrames
using CSV
using Revise
using FaADE


plot = true


θ = 0.5

# Domain
𝒟x = [-0.5, 0.5]
𝒟y = [-0.5, 0.5]


gridoptions = Dict("xbound" => [-0.5, 0.5], "ybound" => [-0.5, 0.5], "xmode" => :stop, "ymode" => :stop)
interpoptions = Dict("interpolant" => :chs)

# Magnetic field
Ψ(x, y) = cos(π * x) * cos(π * y)
# Initial condition
u₀(x, y) = 0.0
# Source term
F(X, t) = 2π^2 * cos(π * X[1]) * cos(π * X[2])


function B(X, x, p, t)
    bn = π * sqrt(abs(cos(x[1] * π) * sin(x[2] * π))^2 + abs(sin(x[1] * π) * cos(x[2] * π))^2)
    X[1] = π * cos(π * x[1]) * sin(π * x[2]) / bn
    X[2] = -π * sin(π * x[1]) * cos(π * x[2]) / bn
    if (x[1] == 0.5) && (x[2] == 0.5)
        X[1] = 0.0
        X[2] = -1.0
    elseif (x[1] == 0.5) && (x[2] == -0.5)
        X[1] = -1.0
        X[2] = 0.0
    elseif (x[1] == -0.5) && (x[2] == -0.5)
        X[1] = 0.0
        X[2] = 1.0
    elseif (x[1] == -0.5) && (x[2] == 0.5)
        X[1] = 1.0
        X[2] = 0.0
    elseif (x[1] == 0.0) && (x[2] == 0.0)
        X[1] = 0.0
        X[2] = 0.0
    end
    # X[3] = 0.0
end
MagField(X, t) = [
    π * cos(π * X[1]) * sin(π * X[2]),
    -π * sin(π * X[1]) * cos(π * X[2]),
    0.0
]

# Exact solution
k_perp = 1.0
T(X, t) = (1.0 - exp(-2.0 * π^2 * k_perp * t)) * Ψ(X[1], X[2]) / k_perp # k_perp = 1


Nx = collect(11:10:41)
Ny = collect(21:20:81)
Exponents = [3.0, 5.0, 6.0, 7.0, 9.0]

coord = :Cartesian

for order in [2, 4]
    dictout = Dict{String,Any}()
    dictout["N"] = Ny

    for EXP in Exponents # Exponent loop
        @show EXP
        k = 1.0 #perpendicular diffusion
        k_para = 10^EXP

        rel_error = []

        for I in eachindex(Nx) # Grid loop
            @show EXP, I, Nx[I]
            nx = Nx[I]
            ny = Ny[I]

            D1 = Grid2D([-0.5, 0.0], [-0.5, 0.5], nx, ny)
            D2 = Grid2D([0.0, 0.5], [-0.5, 0.5], nx, ny)

            joints = ((Joint(2, Right),), (Joint(1, Left),))

            Dom = GridMultiBlock((D1, D2), joints)


            Boundary1Left = SAT_Dirichlet((y, t) -> 0.0, D1.Δx, Left, order)
            Boundary2Right = SAT_Dirichlet((y, t) -> 0.0, D2.Δx, Right, order)

            Boundary1Up = SAT_Dirichlet((x, t) -> 0.0, D1.Δy, Up, order)
            Boundary1Down = SAT_Dirichlet((x, t) -> 0.0, D1.Δy, Down, order)

            Boundary2Up = SAT_Dirichlet((x, t) -> 0.0, D2.Δy, Up, order)
            Boundary2Down = SAT_Dirichlet((x, t) -> 0.0, D2.Δy, Down, order)

            BC = Dict(1 => (Boundary1Left, Boundary1Up, Boundary1Down),
                2 => (Boundary2Right, Boundary2Up, Boundary2Down))




            gdata = construct_grid(B, Dom, [-1.0, 1.0], gridoptions=gridoptions)
            PData = ParallelMultiBlock(gdata, Dom, order, κ=k_para, interpopts=interpoptions)

            P = Problem2D(order, u₀, k, k, Dom, BC, parallel=PData, source=F)

            # Time setup
            Δt = 0.1D1.Δx^2
            t_f = 0.1
            nf = round(t_f / Δt)
            Δt = t_f / nf


            soln = solve(P, Dom, Δt, 2.1Δt)
            soln = solve(P, Dom, Δt, t_f)


            T_exact = [zeros(eltype(Dom.Grids[I]), size(Dom.Grids[I])) for I in eachindex(Dom.Grids)]
            for I in eachindex(Dom.Grids)
                for J in eachindex(Dom.Grids[I])
                    T_exact[I][J] = T(Dom.Grids[I][J], soln.t[2])
                end
            end

            tmp1 = T_exact[1] .- soln.u[2][1]
            tmp2 = T_exact[2] .- soln.u[2][2]

            H = [FaADE.Derivatives.innerH(Dom.Grids[I].Δx, Dom.Grids[I].Δy, Dom.Grids[I].nx, Dom.Grids[I].ny, order) for I in 1:2]
            RE = sqrt(H[1](tmp1, tmp1) + H[2](tmp2, tmp2)) / sqrt(H[1](T_exact[1], T_exact[1]) + H[2](T_exact[2], T_exact[2]))

            push!(rel_error, RE)


        end

        dictout[string("rel ", EXP)] = rel_error

    end

    df = DataFrame(dictout)

    CSV.write("Paper2/data/Convergence/NIMROD_MV_MMS/NIMROD_MultiBlock_O$(order).csv", df)

end
