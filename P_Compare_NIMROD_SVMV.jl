#=
    Compare the NIMROD single block and multiblock cases.
=#

using LinearAlgebra
using Revise
using FaADE


plot = true


θ = 0.5

order = 2

k_perp = 1.0
k_para = 1e9

# Domain
𝒟x = [-0.5, 0.5]
𝒟y = [-0.5, 0.5]


# interpmode = :bilinear
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
T(X, t) = (1.0 - exp(-2.0 * π^2 * k_perp * t)) * Ψ(X[1], X[2]) / k_perp # k_perp = 1

# maxthis = 0.5870063248344728
# texact(t) = (1.0 - exp(-2.0*π^2*t) ) - maxthis
# find_zero(texact,(0.02,0.1),Bisection())

coord = :Cartesian


#=== MIDLINE ===#

nx_off = 11

nx = 11
ny = 41
D11 = Grid2D([-0.5, -0.25], [-0.5, 0.5], nx, ny)
D21 = Grid2D([-0.25, 0.5], [-0.5, 0.5], ny - nx + 1, ny)
DomMV1 = GridMultiBlock((D11, D21), ((Joint(2, Right),), (Joint(1, Left),)))

Boundary11Left = SAT_Dirichlet((y, t) -> 0.0, D11.Δx, Left, order)
Boundary21Right = SAT_Dirichlet((y, t) -> 0.0, D21.Δx, Right, order)

Boundary11Up = SAT_Dirichlet((x, t) -> 0.0, D11.Δy, Up, order)
Boundary11Down = SAT_Dirichlet((x, t) -> 0.0, D11.Δy, Down, order)

Boundary21Up = SAT_Dirichlet((x, t) -> 0.0, D21.Δy, Up, order)
Boundary21Down = SAT_Dirichlet((x, t) -> 0.0, D21.Δy, Down, order)

BCMV1 = Dict(1 => (Boundary11Left, Boundary11Up, Boundary11Down),
    2 => (Boundary21Right, Boundary21Up, Boundary21Down))

gdataMV1 = construct_grid(B, DomMV1, [-1.0, 1.0], gridoptions=gridoptions)
PDataMV1 = ParallelMultiBlock(gdataMV1, DomMV1, order, κ=k_para, interpopts=interpoptions)

PMV1 = Problem2D(order, u₀, k_perp, k_perp, DomMV1, BCMV1, source=F, parallel=PDataMV1)

# Time setup
Δt = 0.1D11.Δx^2
t_f = 0.1
nf = round(t_f / Δt)
Δt = t_f / nf


solnMV1 = solve(PMV1, DomMV1, Δt, 2.1Δt)
solnMV1 = solve(PMV1, DomMV1, Δt, t_f)

T_exact_MV1 = [zeros(size(D11)), zeros(size(D21))]
for I in eachindex(D11)
    T_exact_MV1[1][I] = T(D11[I], t_f)
end
for I in eachindex(D21)
    T_exact_MV1[2][I] = T(D21[I], t_f)
end



#=== OFFSET ===#

nx2_off = 21

nx = 21
ny = 41
D12 = Grid2D([-0.5, 0.0], [-0.5, 0.5], nx, ny)
D22 = Grid2D([0.0, 0.5], [-0.5, 0.5], ny - nx + 1, ny)
DomMV2 = GridMultiBlock((D12, D22), ((Joint(2, Right),), (Joint(1, Left),)))

Boundary12Left = SAT_Dirichlet((y, t) -> 0.0, D12.Δx, Left, order)
Boundary22Right = SAT_Dirichlet((y, t) -> 0.0, D22.Δx, Right, order)

Boundary12Up = SAT_Dirichlet((x, t) -> 0.0, D12.Δy, Up, order)
Boundary12Down = SAT_Dirichlet((x, t) -> 0.0, D12.Δy, Down, order)

Boundary22Up = SAT_Dirichlet((x, t) -> 0.0, D22.Δy, Up, order)
Boundary22Down = SAT_Dirichlet((x, t) -> 0.0, D22.Δy, Down, order)

BCMV2 = Dict(1 => (Boundary12Left, Boundary12Up, Boundary12Down),
    2 => (Boundary22Right, Boundary22Up, Boundary22Down))

gdataMV2 = construct_grid(B, DomMV2, [-1.0, 1.0], gridoptions=gridoptions)
PDataMV2 = ParallelMultiBlock(gdataMV2, DomMV2, order, κ=k_para, interpopts=interpoptions)

PMV2 = Problem2D(order, u₀, k_perp, k_perp, DomMV2, BCMV2, source=F, parallel=PDataMV2)

# Time setup
Δt = 0.1D12.Δx^2
t_f = 0.1
nf = round(t_f / Δt)
Δt = t_f / nf


solnMV2 = solve(PMV2, DomMV2, Δt, 2.1Δt)
solnMV2 = solve(PMV2, DomMV2, Δt, t_f)

T_exact_MV2 = [zeros(size(D12)), zeros(size(D22))]
for I in eachindex(D12)
    T_exact_MV2[1][I] = T(D12[I], t_f)
end
for I in eachindex(D22)
    T_exact_MV2[2][I] = T(D22[I], t_f)
end







# Single volume run for reference solution

nx = 41
ny = 41
Dom = Grid2D(𝒟x, 𝒟y, nx, ny)

BoundaryLeft = SAT_Dirichlet((y, t) -> 0.0, Dom.Δx, Left, order)
BoundaryRight = SAT_Dirichlet((y, t) -> 0.0, Dom.Δx, Right, order)
BoundaryUp = SAT_Dirichlet((x, t) -> 0.0, Dom.Δy, Up, order)
BoundaryDown = SAT_Dirichlet((x, t) -> 0.0, Dom.Δy, Down, order)

BC = (BoundaryLeft, BoundaryRight, BoundaryUp, BoundaryDown)

gdata = construct_grid(B, Dom, [-1.0, 1.0], ymode=:stop)
PData = ParallelData(gdata, Dom, order, κ=k_para, interpolant=:chs)

P = Problem2D(order, u₀, k_perp, k_perp, Dom, BC, source=F, parallel=PData)

# Time setup
Δt = 0.1Dom.Δx^2
t_f = 0.1
nf = round(t_f / Δt)
Δt = t_f / nf


soln = solve(P, Dom, Δt, 2.1Δt)
soln = solve(P, Dom, Δt, t_f)

T_exact = zeros(size(Dom))
for I in eachindex(Dom.gridx)
    T_exact[I] = T(Dom[I], t_f)
end





println("plotting")

using GLMakie
# using CairoMakie



# absolute errors
region11 = abs.(solnMV1.u[2][1] .- T_exact_MV1[1]) ./ T_exact_MV1[1];
region21 = abs.(solnMV1.u[2][2] .- T_exact_MV1[2]) ./ T_exact_MV1[2];
region12 = abs.(solnMV2.u[2][1] .- T_exact_MV2[1]) ./ T_exact_MV2[1];
region22 = abs.(solnMV2.u[2][2] .- T_exact_MV2[2]) ./ T_exact_MV2[2];


# Colour range for first plot
tmpcr1 = (minimum(minimum.(region11[2:end, 2:end-1])), maximum(maximum.(region11[2:end, 2:end-1])))
tmpcr2 = (minimum(minimum.(region21[1:end-1, 2:end-1])), maximum(maximum.(region21[1:end-1, 2:end-1])))

colourrange1 = (min(tmpcr1[1], tmpcr2[1]), max(tmpcr1[2], tmpcr2[2]))

# Colour range for second plot
tmpcr1 = (minimum(minimum.(region12[2:end, 2:end-1])), maximum(maximum.(region12[2:end, 2:end-1])))
tmpcr2 = (minimum(minimum.(region22[1:end-1, 2:end-1])), maximum(maximum.(region22[1:end-1, 2:end-1])))

colourrange2 = (min(tmpcr1[1], tmpcr2[1]), max(tmpcr1[2], tmpcr2[2]))






f = Figure(size=(1200, 500), fontsize=20)
axgg = f[1, 1] = GridLayout()

ax1 = Axis(axgg[1, 1], xlabel="x", ylabel="y", xticklabelrotation=pi / 4, xlabelpadding=-5.0)
ax2 = Axis(axgg[1, 3], xlabel="x", ylabel="y", xticklabelrotation=pi / 4, xlabelpadding=-5.0)




###### First plot - grid centred
println("First plot")

gax21 = surface!(ax1, DomMV2.Grids[1].gridx[2:end, 2:end-1], DomMV2.Grids[1].gridy[2:end, 2:end-1], region12[2:end, 2:end-1], colorrange=colourrange2, colormap=:viridis, shading=NoShading)
gax22 = surface!(ax1, DomMV2.Grids[2].gridx[1:end-1, 2:end-1], DomMV2.Grids[2].gridy[1:end-1, 2:end-1], region22[1:end-1, 2:end-1], colorrange=colourrange2, colormap=:viridis, shading=NoShading)

wireframe!(ax1, DomMV2.Grids[1].gridx, DomMV2.Grids[1].gridy, region12, color=(:black, 0.1), overdraw=true)
wireframe!(ax1, DomMV2.Grids[2].gridx, DomMV2.Grids[2].gridy, region22, color=(:black, 0.1), overdraw=true)

lines!(ax1, DomMV2.Grids[1].gridx[end, :], DomMV2.Grids[1].gridy[end, :], region11[end, :], color=:magenta, linestyle=:dash, linewidth=2.0, overdraw=true)


###### Second plot - grid offset
println("Second plot")


gax11 = surface!(ax2, DomMV1.Grids[1].gridx[2:end, 2:end-1], DomMV1.Grids[1].gridy[2:end, 2:end-1], region11[2:end, 2:end-1], colorrange=colourrange1, colormap=:viridis, shading=NoShading)
gax12 = surface!(ax2, DomMV1.Grids[2].gridx[1:end-1, 2:end-1], DomMV1.Grids[2].gridy[1:end-1, 2:end-1], region21[1:end-1, 2:end-1], colorrange=colourrange1, colormap=:viridis, shading=NoShading)

wireframe!(ax2, DomMV1.Grids[1].gridx, DomMV1.Grids[1].gridy, region11, color=(:black, 0.1), overdraw=true)
wireframe!(ax2, DomMV1.Grids[2].gridx, DomMV1.Grids[2].gridy, region21, color=(:black, 0.1), overdraw=true)

lines!(ax2, DomMV1.Grids[1].gridx[end, :], DomMV1.Grids[1].gridy[end, :], region12[end, :], color=:magenta, linestyle=:dash, linewidth=2.0, overdraw=true)





# gaxoo = surface!(axo,Dom.gridx[2:end-1,2:end-1], Dom.gridy[2:end-1,2:end-1], region[2:end-1,2:end-1], colormap=:viridis, shading=NoShading)
# wireframe!(axo, Dom.gridx, Dom.gridy, region, color=(:black,0.1), overdraw=true)


Colorbar(axgg[1, 2], limits=colourrange2, ticklabelsize=20, labelsize=20)
Colorbar(axgg[1, 4], limits=colourrange1, label=L"|T_{\text{MV}} - T_{\text{exact,MV}}|/T_{\text{exact}}", ticklabelsize=20, labelsize=20)
# Colorbar(axgg[1,2],gaxoo, ticklabelsize=20)



xlims!(ax1, (-0.5, 0.5))
ylims!(ax1, (-0.5, 0.5))

xlims!(ax2, (-0.5, 0.5))
ylims!(ax2, (-0.5, 0.5))


hideydecorations!(ax2)
# hideydecorations!(ax1)





println("saving")

# f
# save("2D/Figures/Compare_NIMROD_SVMV.svg", f)
save("Paper2/data/Compare_NIMROD_SVMV_linear.png", f, px_per_unit=4.0)
