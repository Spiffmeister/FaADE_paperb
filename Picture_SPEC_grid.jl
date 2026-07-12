
using FaADE
using GLMakie
using SPECReader





θ = 0.5
order = 2

k_para = 1.0e6
k_perp = 1.0


# Time setup
Δt = 1e-4
t_f = 1e-2
nf = round(t_f / Δt)
Δt = t_f / nf


gridoptions = Dict("xbound" => [-0.5, 0.5], "ybound" => [-0.5, 0.5], "xmode" => :stop, "ymode" => :stop)
interpoptions = Dict("interpolant" => :chs)

# Magnetic field
Ψ(x, y) = cos(π * x) * cos(π * y)
# Initial condition
u₀(x, y) = 0.0


coord = :Curvilinear







speceq = SPECEquilibrium("./data/G3V01L0Fi.002.sp.h5")

specboundary = get_boundary(speceq, 1)
specaxis = get_axis(speceq)

Tor = FaADE.Grid.Torus(specboundary.Rbc, specboundary.Zbs, specboundary.m, specboundary.n)


lrdil = [0.4, 0.6]
tddil = 0.0

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



NGrids = [(11, 11),
    (11, 11),
    (11, 11),
    (11, 11),
    (11, 11)]



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

D2 = Grid2D(d2b, d2l, d2r, d2t, NGrids[2]...)


# Top domain
d3b = d1t
d3l(v) = (TorCentre + [-width_l, height_l]) * (1 - v) + Tor(arcs[2], 0.0) * v
d3r(v) = (TorCentre + [width_r, height_r]) * (1 - v) + Tor(arcs[1], 0.0) * v
d3t(u) = Tor(u * (-arcs[2] + arcs[1]) + arcs[2], 0.0)

D3 = Grid2D(d3b, d3l, d3r, d3t, NGrids[3]...)



# Left domain
d4b = d1l
d4l(v) = (TorCentre + [-width_l, -height_l]) * (1 - v) + v * Tor(arcs[3], 0.0)
d4r(v) = (TorCentre + [-width_l, height_l]) * (1 - v) + v * Tor(arcs[2], 0.0)
d4t(u) = Tor(u * (-arcs[3] + arcs[2]) + arcs[3], 0.0)

D4 = Grid2D(d4b, d4l, d4r, d4t, NGrids[4]...)



# Bottom domain - bottom is right-left of d1b
d5b(u) = TorCentre + [-width_l, -height_l] + u * ([width_r, -height_r] - [-width_l, -height_l]) + u * (1 - u) * [0.0, -tddil]
d5l(v) = (TorCentre + [width_r, -height_r]) * (1 - v) + v * Tor(arcs[4], 0.0)
d5r(v) = (TorCentre + [-width_l, -height_l]) * (1 - v) + v * Tor(arcs[3], 0.0)
d5t(u) = Tor(u * (-arcs[4] + arcs[3]) + arcs[4], 0.0)

D5 = Grid2D(d5b, d5l, d5r, d5t, NGrids[5]...)


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






dH(X, x, params, t) = field_line!(X, t, x, speceq)

XtoB(x, y) = find_sθζ((x, y), 0.0, speceq, 1)
BtoX(r, θ) = get_RZ(r, θ, 0.0, speceq, 1)

gridoptions = Dict("coords" => (XtoB, BtoX))

# gdata = construct_grid(dH, Dom, [-2π, 2π], gridoptions=gridoptions)





# Source term
inset = 0.90
function source(X, t)
    s, _ = find_sθζ(X, 0.0, speceq, 1) # This is cubersome and should be replaced, maybe by a nearest neighbour interpolation
    s = (s + 1) / 2 # s ∈ [0,1]
    tmp = exp(-s^2 / 0.4^2)
    if (inset - s) > 0
        return tmp
    else
        return 0
    end
end



f = Figure(fontsize=20, size=(500, 600))
axf = Axis(f[1, 1], xlabel="R", ylabel="Z")
for I in 1:5
    wireframe!(axf, Dom.Grids[I].gridx, Dom.Grids[I].gridy, zeros(size(Dom.Grids[I].gridx)), linewidth=1)

    lines!(axf, Dom.Grids[I].gridx[1, :], Dom.Grids[I].gridy[1, :], color=(:red), overdraw=true)
    lines!(axf, Dom.Grids[I].gridx[end, :], Dom.Grids[I].gridy[end, :], color=(:red), overdraw=true)
    lines!(axf, Dom.Grids[I].gridx[:, 1], Dom.Grids[I].gridy[:, 1], color=(:red), overdraw=true)
    lines!(axf, Dom.Grids[I].gridx[:, end], Dom.Grids[I].gridy[:, end], color=(:red), overdraw=true)

end


poinout = ReadPoincare("./data/G3V01L0Fi.002.sp.h5")

scatter!(axf, poinout.R[1, :, :][:], poinout.Z[1, :, :][:], markersize=2, color=(:black))

ylims!(axf, minimum(Dom.Grids[5].gridy), maximum(Dom.Grids[3].gridy))
xlims!(axf, minimum(Dom.Grids[3].gridx), maximum(Dom.Grids[2].gridx))


save("./data/SPEC_GridPoincare.png", f, px_per_unit=4)


f
