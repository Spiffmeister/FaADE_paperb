#=
    Picture of the single island error in a circle
=#
using LinearAlgebra
using Revise
using CubicHermiteSpline
using Surrogates
using FaADE


using JLD2
using DataFrames
using CSV



save_plots = true


### Fixed parameters

order = 2
EXP = 6
K_para = 10.0^EXP
K_perp = 1.0

t_f = 1e-2

resolutions = 21:10:81
reference_resolution = 101

# reference_dt = 1e-5
reference_dt = nothing


dilation = 0.0 #dilation of the interior domain
g0b = 0.25
α = 0.1

if α == 0.0
    packfn(x) = x
else
    island_location = (0.7 - (g0b + 0.25 * dilation)) / (1.0 - (g0b + 0.25 * dilation))
    Dp(x, α, B, B0, B1) = B + α * sinh(asinh((B1 - B) / α) * x + asinh((B0 - B) / α) * (1 - x))
    packfn(x) = Dp(x, α, island_location, 0.0, 1.0)
end



filepaths = [
    string("Paper2/data/Convergence/SingleIslandSelf_2_6_0.1_0.25_0.1"),
    string("Paper2/data/Convergence/SingleIslandSelf_2_9_0.1_0.25_0.1")
]



EXP = string.("\\kappa_\\parallel=", ["10^6", "10^9"])



println("Start")



N = [(n, n) for n in resolutions]
Doms = [FaADE.Grid.squared_circle([n for _ in 1:5], radial_packing=packfn) for n in N]



nx = ny = reference_resolution
reference_Dom = FaADE.Grid.squared_circle([(nx, ny) for _ in 1:5], radial_packing=packfn)
# @load string(dirname,"/ref") reference_u_O2




# using GLMakie
using CairoMakie

# Solution across the separatrix
f = Figure(fontsize=16)
axf = [Axis(f[i, 1], xlabel=L"x", ylabel=L"u(x,0,t_f)") for i in eachindex(filepaths)]

for (fi, filepath) in enumerate(filepaths)
    file = jldopen(string(filepath, "/runs"))
    solns_u = file["solns_u"]

    # reference_Dom = FaADE.Grid.squared_circle([(nx,ny) for _ in 1:5], radial_packing=packfn)
    @load string(filepath, "/ref") reference_u


    # for I in 1:3:length(solns_u)
    li = [
        lines!(axf[fi],
            Doms[I].Grids[4].gridx[floor(Int, Doms[I].Grids[4].nx / 2 + 1), :],
            solns_u[I][4][floor(Int, Doms[I].Grids[4].nx / 2 + 1), :],
            label=string("n=", Doms[I].Grids[4].nx)) for I in 1:3:length(solns_u)
    ]
    # end

    xlims!(-0.9, -0.5)
    ylims!(0.0, 0.004)

    text!(axf[fi], 0, 1, text=L"%$(EXP[fi])", space=:relative, align=(:left, :top), offset=(1, 0), fontsize=22)
    # axislegend(axf[fi])
    if fi == 1
        Legend(f[1, 2], li, string.([Doms[I].Grids[4].nx for I in 1:3]))
    end
end

hidexdecorations!.(axf[1], grid=false)
linkaxes!(axf...)


f

if save_plots
    save("lines.pdf", f)
end



# ==== TRICONTOURF ==== #


δ = 0.05
rs = 0.7
function B(X, x::Array{Float64}, params, t)
    X[1] = δ * x[1] * (1 - x[1]) * sin(x[2])#/bn
    X[2] = (2x[1] - 2 * rs + δ * (1 - x[1]) * cos(x[2]) - δ * x[1] * cos(x[2]))#/bn
end
dH(X, x, params, t) = B(X, x, params, t)

include("../FieldLines.jl")
poindata = FieldLines.construct_poincare(dH, [0.0, 1.0], [0.0, π], N_trajs=400, N_orbs=400)
BtoX(r, θ) = [r * cos(θ), r * sin(θ)]
poinrtheta = hcat([BtoX(poindata.ψ[I], poindata.θ[I]) for I in eachindex(poindata.ψ)]...)



g = Figure(size=(8, 7) .* 120, fontsize=16)
# axg = Axis(g[1,1], xlabel="x", ylabel="y", aspect=DataAspect())
axg = Axis(g[1, 1], xlabel="x", ylabel="y")

file = jldopen(string(filepaths[2], "/ref"))
reference_u = file["reference_u"]

gridx = vcat([reference_Dom.Grids[I].gridx for I in 1:5]...)[:]
gridy = vcat([reference_Dom.Grids[I].gridy for I in 1:5]...)[:]
utot = vcat(reference_u...)[:]

utot = utot[gridy.≥-1e-11]
gridx = gridx[gridy.≥-1e-11]
gridy = gridy[gridy.≥-1e-11]

tcf = tricontourf!(axg, gridx, gridy, utot, levels=0.0:0.002/10:0.002)


# file = jldopen(string(filepaths[2],"/ref"))
# reference_u = file["reference_u"]

# gridx = vcat([reference_Dom.Grids[I].gridx for I in 1:5]...)[:]
# gridy = vcat([reference_Dom.Grids[I].gridy for I in 1:5]...)[:]
# utot = vcat(reference_u...)[:]

# utot = utot[gridy .≤ 1e-11]
# gridx = gridx[gridy .≤ 1e-11]
# gridy = gridy[gridy .≤ 1e-11]

# tcf = tricontourf!(axg, gridx, gridy, utot, levels=0.0:0.002/10:0.002)

Colorbar(g[1, 2], tcf, height=Relative(0.9), label=L"u(x,y,t_f)", labelsize=18)

xlims!(-1.0, 1.0)
ylims!(-0.0, 1.0)



poinrmask = 0.55 .< poindata.ψ .< 0.85
scatter!(axg, poinrtheta[1, poinrmask], poinrtheta[2, poinrmask], markersize=1.5, color=(:red, 0.9))

hidexdecorations!(axg, grid=false)

g

if save_plots
    save(string(filepaths[2], "/island_tricontour.png"), g, px_per_unit=4)
end





# ==== Plot local relative error ==== #



function reconstruct_soln(Dom, refinterp)
    soln = [zeros(Dom.Grids[I].nx, Dom.Grids[I].ny) for I in 1:5]
    for I in eachindex(Dom.Grids)
        for J in eachindex(Dom.Grids[I].gridx)
            # soln[I][J] = evaluate(refinterp[I],[Dom.Grids[I][J]...])[1]
            soln[I][J] = refinterp[I](tuple(Dom.Grids[I][J]...))
        end
    end
    return soln
end

function compute_errors(ublock, refu, Dom)
    tmpdenom = 0.0
    for I in eachindex(ublock) #for each grid

        H = FaADE.Derivatives.innerH(Dom.Grids[I].Δx, Dom.Grids[I].Δy, Dom.Grids[I].nx, Dom.Grids[I].ny, order)

        tmpdenom += H(refu[I], Dom.Grids[I].J, refu[I])
    end
    return sqrt(tmpdenom)
end



refinterp = []
for I in eachindex(reference_Dom.Grids)
    grid = tuple.(reference_Dom.Grids[I].gridx[:], reference_Dom.Grids[I].gridy[:])
    push!(refinterp, RadialBasis(grid, reference_u[I], minimum(reference_u[I]), maximum(reference_u[I]), rad=cubicRadial()))
end

file = jldopen(string(filepaths[1], "/runs"))
solns_u = file["solns_u"]

soln_ref = reconstruct_soln(Doms[5], refinterp)

# soln_err = [abs.(solns) for solns in soln_ref .- solns_u[5]]/abs(maximum(maximum.([abs.(solns) for solns in soln_ref .- solns_u[5]])))

soln_err = [abs.(solns) for solns in soln_ref .- solns_u[5]] / compute_errors(solns_u[5], soln_ref, Doms[5])


# h = Figure(fontsize=18, size=(10,5).*120)
# axh = Axis(h[1,1], xlabel="x", ylabel="y", aspect=DataAspect())

# axh = Axis(g[2,1], xlabel="x", ylabel="y", aspect=DataAspect())
axh = Axis(g[2, 1], xlabel="x", ylabel="y")

cmap1 = minimum(minimum.(soln_err)), maximum(maximum.(soln_err))


for I in 1:5
    surface!(axh, Doms[5].Grids[I].gridx, Doms[5].Grids[I].gridy, soln_err[I], colorrange=cmap1, colorscale=log)
end
# Colorbar(h[1,2], limits=cmap1, label=L"\log\left((u_{ref} - u)/\max(u_{ref} - u)\right)",height=Relative(0.95))
Colorbar(g[2, 2], limits=cmap1, label=L"\log\left((u_{ref} - u)/||u_{ref}||_H\right)", height=Relative(0.95))

xlims!(-1.0, 1.0)
ylims!(-1.0, 0.0)


colgap!(g.layout, 10)


hidexdecorations!(axh, grid=false)



if save_plots
    save(string(filepaths[1], "/island_relerror.png"), g, px_per_unit=4)
end
