
using FaADE
using LinearAlgebra
using JLD2
using DelaunayTriangulation
using CubicHermiteSpline
using SPECReader

using GLMakie
using CairoMakie

GLMakie.activate!()




compute_run = false
compute_pmap = false
save_run = true
save_pmap = true
save_figure = true


order = 2

k_para = 1.0e9
k_perp = 1.0e0


# Time setup
Δt = 1e-8
t_f = 5e-3


# Initial condition
u₀(x, y) = 0.0


coord = :Curvilinear






speceq = SPECEquilibrium("data/G3V01L0Fi.002.sp.h5")

specboundary = get_boundary(speceq, 1)
specaxis = get_axis(speceq)

@show "CREATE GRIDS"


Tor = FaADE.Grid.Torus(specboundary.Rbc, specboundary.Zbs, specboundary.m, specboundary.n)


lrdil = [0.4, 0.6]
tddil = 0.2

shift = 1.e-1
width = 0.1
width_l = 0.075
width_r = 0.075
height_l = 0.35
height_r = 0.5
TorCentre = [sum(specaxis.R[1:2]) - shift, 0.0]


arcs = [π / 4 + π / 24, 3π / 4 + π / 12, 5π / 4 - π / 12, 7π / 4 - π / 24]
arcs = vcat(arcs, arcs[1] + 2π)
arcs = -arcs



NGrids = [(61, 61),
    (61, 71),
    (61, 71),
    (61, 71),
    (61, 71)]
# NGrids = [(101, 101) for _ in 1:5]

Dp(x, α, B, B0, B1) = B + α * sinh(asinh((B1 - B) / α) * x + asinh((B0 - B) / α) * (1 - x))
packfn(x) = Dp(x, 0.3, 0.65, 0.0, 1.0)

# Mid domain
d1b(u) = TorCentre + [-width_l, -height_l] + u * ([width_r, -height_r] - [-width_l, -height_l]) + u * (1 - u) * [0.0, -tddil]
d1l(v) = TorCentre + [-width_l, -height_l] + v * ([-width_l, height_l] - [-width_l, -height_l]) + v * (1 - v) * [lrdil[1], 0.0]
d1r(v) = TorCentre + [width_r, -height_r] + v * ([width_r, height_r] - [width_r, -height_r]) + v * (1 - v) * [lrdil[2], 0.0]
d1t(u) = TorCentre + [-width_l, height_l] + u * ([width_r, height_r] - [-width_l, height_l]) + u * (1 - u) * [0.0, tddil]

D1 = Grid2D(d1b, d1l, d1r, d1t, NGrids[1]...)



# Right domain
d2b(u) = TorCentre + [width_r, height_r] + u * ([width_r, -height_r] - [width_r, height_r]) + u * (1 - u) * [lrdil[2], 0.0]
d2l(v) = (TorCentre + [width_r, height_r]) * (1 - v) + Tor(arcs[end], 0.0) * v
d2r(v) = (TorCentre + [width_r, -height_r]) * (1 - v) + Tor(arcs[end-1], 0.0) * v
d2t(u) = Tor(u * (-arcs[end] + arcs[end-1]) + arcs[end], 0.0)

D2 = Grid2D(d2b, d2l, d2r, d2t, NGrids[2]..., stretchv=packfn)


# Top domain
d3b = d1t
d3l(v) = (TorCentre + [-width_l, height_l]) * (1 - v) + Tor(arcs[2], 0.0) * v
d3r(v) = (TorCentre + [width_r, height_r]) * (1 - v) + Tor(arcs[1], 0.0) * v
d3t(u) = Tor(u * (-arcs[2] + arcs[1]) + arcs[2], 0.0)

D3 = Grid2D(d3b, d3l, d3r, d3t, NGrids[3]..., stretchv=packfn)



# Left domain
d4b = d1l
d4l(v) = (TorCentre + [-width_l, -height_l]) * (1 - v) + v * Tor(arcs[3], 0.0)
d4r(v) = (TorCentre + [-width_l, height_l]) * (1 - v) + v * Tor(arcs[2], 0.0)
d4t(u) = Tor(u * (-arcs[3] + arcs[2]) + arcs[3], 0.0)

D4 = Grid2D(d4b, d4l, d4r, d4t, NGrids[4]..., stretchv=packfn)



# Bottom domain - bottom is right-left of d1b
d5b(u) = TorCentre + [-width_l, -height_l] + u * ([width_r, -height_r] - [-width_l, -height_l]) + u * (1 - u) * [0.0, -tddil]
d5l(v) = (TorCentre + [width_r, -height_r]) * (1 - v) + v * Tor(arcs[4], 0.0)
d5r(v) = (TorCentre + [-width_l, -height_l]) * (1 - v) + v * Tor(arcs[3], 0.0)
d5t(u) = Tor(u * (-arcs[4] + arcs[3]) + arcs[4], 0.0)

D5 = Grid2D(d5b, d5l, d5r, d5t, NGrids[5]..., stretchv=packfn)


joints = ((Joint(2, FaADE.Right), Joint(3, FaADE.Up), Joint(4, FaADE.Left), Joint(5, FaADE.Down)),
    (Joint(1, FaADE.Down), Joint(3, FaADE.Left), Joint(5, FaADE.Right)),
    (Joint(1, FaADE.Down), Joint(4, FaADE.Left), Joint(2, FaADE.Right)),
    (Joint(1, FaADE.Down), Joint(5, FaADE.Left), Joint(3, FaADE.Right)),
    (Joint(1, FaADE.Down), Joint(2, FaADE.Left), Joint(4, FaADE.Right)))

Dom = GridMultiBlock((D1, D2, D3, D4, D5), joints)





Bxy(X, t) = 0.0

# Boundary conditions
Dr = SAT_Dirichlet(Bxy, D2.Δy, FaADE.Up, order, D2.Δx, :Curvilinear) # Block 2 BCs
Du = SAT_Dirichlet(Bxy, D3.Δy, FaADE.Up, order, D3.Δx, :Curvilinear) # Block 3 BCs
Dl = SAT_Dirichlet(Bxy, D4.Δy, FaADE.Up, order, D4.Δx, :Curvilinear) # Block 4 BCs
Dd = SAT_Dirichlet(Bxy, D5.Δy, FaADE.Up, order, D5.Δx, :Curvilinear) # Block 5 BCs

BD = Dict(2 => (Dr,), 3 => (Du,), 4 => (Dl,), 5 => (Dd,))



@show "PARALLEL MAP"


dH(X, x, params, t) = field_line!(X, t, x, speceq)

XtoB(x, y) = find_sθζ((x, y), 0.0, speceq, 1)
BtoX(r, θ) = get_RZ(r, θ, 0.0, speceq, 1)

gridoptions = Dict("coords" => (XtoB, BtoX))



intercept(u, x, y, t) = begin
    if isnan(u)
        return zero(typeof(x))
    else
        return u
    end
end

savedir_grid = "./data/SPEC_Case_n$(D1.nx)_$(D2.ny)/"
savedir = "./data/SPEC_Case_n$(D1.nx)_$(D2.ny)_$(Δt)_$(t_f)/"
if !isdir(savedir)
    mkdir(savedir)
end
if !isdir(savedir_grid)
    mkdir(savedir_grid)
end
if compute_pmap
    println("Computing parallel map")
    gdata = construct_grid(dH, Dom, [-2π, 2π], gridoptions=gridoptions)
    if save_pmap
        for I in eachindex(gdata)
            jldsave(string(savedir_grid, "pgrid_", I);
                pgrid=gdata[I]
            )
        end
    end
else
    println("Reading parallel map")
    gdata = Dict()
    for I in 1:5
        pgrid_file = jldopen(string(savedir_grid, "pgrid_", I))
        gdata[I] = pgrid_file["pgrid"]
    end
end


function magfield(X, t)
    s, θ = find_sθζ(X, 0.0, speceq, 1)
    return get_Bfield(s, θ, 0.0, speceq)
end

interpotions = Dict("interpolant" => :chs, "intercept" => intercept)


PData = FaADE.ParallelOperator.ParallelMultiBlock(gdata, Dom, order, κ=k_para, interpopts=interpotions)


@show "SOURCE CONSTRUCTION"



inset = 0.9
function source_interp(X)
    s = find_sθζ(X, 0.0, speceq, 1)[1]
    s = (s + 1) / 2 # s ∈ [-1,1] → [0,1]
    tmp = 4 * (1 - s^2)^8
    if (inset - s) > 0
        return tmp
    else
        return zero(eltype(X))
    end
end



source = [[source_interp(grid[i, j]) for i in 1:grid.nx, j in 1:grid.ny] for grid in Dom.Grids]



@show "SOLVE"

# Build PDE problem
P = Problem2D(order, u₀, k_perp, k_perp, Dom, BD, source=source, parallel=PData)

nf = round(t_f / Δt)
Δt = t_f / nf

if compute_run
    solve(P, Dom, Δt, Δt)
    soln = solve(P, Dom, Δt, t_f)

    soln_u = soln.u[2]

    if save_run
        jldsave(string(savedir, "run");
            soln_u=soln.u[2]
        )
    end
else
    file = jldopen(string(savedir, "run"))
    soln_u = file["soln_u"]
end


@show "SOLVE DONE"


#=== PLOTS ===#

cmap = minimum(minimum.(soln_u)), maximum(maximum.(soln_u))

# Plot grid and Poincare
f = Figure()
axf = Axis(f[1, 1])
axf2 = Axis(f[1, 2])
for I in 1:5
    wireframe!(axf, Dom.Grids[I].gridx, Dom.Grids[I].gridy, zeros(size(Dom.Grids[I].gridx)), alpha=0.2)

    surface!(axf2, Dom.Grids[I].gridx, Dom.Grids[I].gridy, soln_u[I], colorrange=cmap)
    wireframe!(axf2, Dom.Grids[I].gridx, Dom.Grids[I].gridy, soln_u[I])
end


poinout = ReadPoincare("./data/G3V01L0Fi.002.sp.h5")

scatter!(axf, poinout["R"][1, :, :][:], poinout["Z"][1, :, :][:], markersize=4, color=(:red), overdraw=true)



# Surface plot
g = Figure()
axg = Axis3(g[1, 1])
for I in 1:5
    surface!(axg, Dom.Grids[I].gridx, Dom.Grids[I].gridy, soln_u[I])
    wireframe!(axg, Dom.Grids[I].gridx, Dom.Grids[I].gridy, soln_u[I])
end

if save_figure
    GLMakie.activate!()
    save(string(savedir, "grid_surface.png"), f, px_per_unit=4)
end





tmpf = Figure();
tmpfax = Axis3(tmpf[1, 1]);
for I in 1:5
    surface!(tmpfax, Dom.Grids[I].gridx, Dom.Grids[I].gridy, source[I])
    wireframe!(tmpfax, Dom.Grids[I].gridx, Dom.Grids[I].gridy, source[I])
end




# clevels = [0.0, maximum(maximum.(soln_u))]
clevels = [2e-4, 7e-4]
nlevels = 11

cmap = clevels[1]:diff(clevels)[1]/nlevels:clevels[2]

contf = Figure(size=(10, 4) .* 150, fontsize=16);
contax = Axis(contf[1, 1], xlabel="R", ylabel="Z");
contour_plots = Contour[]
for I in 1:5
    contplt = contour!(
        contax,
        Dom.Grids[I].gridx,
        Dom.Grids[I].gridy,
        soln_u[I],
        levels=cmap,
        colorrange=(clevels[1], clevels[2]),
        linewidth=2.0
    )
    push!(contour_plots, contplt)
end
scatter!(contax, poinout["R"][1, :, :][:], poinout["Z"][1, :, :][:], markersize=4, color=(:red), overdraw=true)
ylims!(0.0, 0.85)
xlims!(5.525, 6.215)

Colorbar(contf[1, 2],
    limits=clevels,
    colormap=cgrad(:viridis, length(cmap), categorical=true),
    label=L"u(R,Z,t_f)"
)

# contf

if save_figure
    CairoMakie.activate!()
    save(string(savedir, "SPECContour.png"), contf, px_per_unit=4)
end
