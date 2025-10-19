#=
    Pollution in NIMROD becnhmark
=#
using LinearAlgebra
using DelimitedFiles
using CSV
using DataFrames
using Statistics
using BasicInterpolators

using FaADE



θ = 0.5

Ψ(x, y) = cos(π * x) * cos(π * y)


gridoptions = Dict("xbound" => [-0.5, 0.5], "ybound" => [-0.5, 0.5], "xmode" => :stop, "ymode" => :stop)
interpoptions = Dict("interpolant" => :chs)


# Initial condition
u₀(x, y) = 0.0
# Source term
F(X, t) = 2π^2 * cos(π * X[1]) * cos(π * X[2])
# Magnetic field for FLT
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


M = [0.1 / 2^i for i in 0:1]
KPEXP = [0, 3, 5, 6, 7, 9, 10]
KP = 10.0 .^ KPEXP

nx = 21 #fixed grid resolution
ny = 41

coord = :Cartesian

for order in [2, 4]
    k_para = 1.0
    k = k_perp = 1.0
    T(x, y, t) = (1.0 - exp(-2.0 * k_perp * π^2 * t)) / (k_perp) * Ψ(x, y)
    # Diffusion coefficient

    dictout = Dict{String,Any}()
    dictout["KP"] = KPEXP

    for dt in M
        pollution = []
        rel_error = []
        abs_error = []

        comp_error = []
        comp_poll = []

        println(" --- order=", order, " --- θ=", θ, " ---")
        for k_para in KP

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

            # Time setup
            Δt = dt
            t_f = 1.0
            nf = round(t_f / Δt)
            Δt = t_f / nf

            gdata = construct_grid(B, Dom, [-1.0, 1.0], gridoptions=gridoptions)
            PData = ParallelMultiBlock(gdata, Dom, order, κ=k_para, interpopts=interpoptions)

            # Build PDE problem
            P = Problem2D(order, u₀, k, k, Dom, BC, source=F, parallel=PData)

            soln = solve(P, Dom, Δt, 1.1Δt)
            soln = solve(P, Dom, Δt, t_f)

            # Solution without parallel operator
            Pwo = Problem2D(order, u₀, k, k, Dom, BC, source=F, parallel=nothing)
            solnwo = solve(P, Dom, Δt, t_f)

            # T_exact = zeros(eltype(Dom),size(Dom));
            # T_exact = [zeros(size(Dom.Grids[1])), zeros(size(Dom.Grids[2]))]
            # for I in eachindex(Dom)
            #     T_exact[I] = T(Dom.Grids[I]...,soln.t[end])
            # end

            push!(pollution, abs(1 / soln.u[2][1][end, floor(Int, ny / 2)+1] - 1))

        end

        dictout[string("poll ", dt)] = pollution

    end
    df = DataFrame(dictout)
    CSV.write(string("Paper2/data/Convergence/NIMROD_MV_MMS/NBMVAP_O$(order).csv"), df)
end
