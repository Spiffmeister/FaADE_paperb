using LinearAlgebra

using FaADE


n = 11
dilation = 0.1


Tor = FaADE.Grid.Torus([1.0], [1.0], [1], [0])

D1 = Grid2D(
    u -> [-0.25, -0.25] + u * ([0.25, -0.25] - [-0.25, -0.25]) + u * (1 - u) * [0.0, -dilation],
    v -> [-0.25, -0.25] + v * ([-0.25, 0.25] - [-0.25, -0.25]) + v * (1 - v) * [-dilation, 0.0],
    v -> [0.25, -0.25] + v * ([0.25, 0.25] - [0.25, -0.25]) + v * (1 - v) * [dilation, 0.0],
    u -> [-0.25, 0.25] + u * ([0.25, 0.25] - [-0.25, 0.25]) + u * (1 - u) * [0.0, dilation],
    n, n)

D2 = Grid2D(#u->[0.25, -u*0.5 + 0.25], # Bottom
    u -> [0.25, 0.25] + u * ([0.25, -0.25] - [0.25, 0.25]) + u * (1 - u) * [dilation, 0.0],
    v -> v * (Tor(π / 4, 0.0) - [0.25, 0.25]) + [0.25, 0.25], # Left
    v -> v * (Tor(7π / 4, 0.0) + [-0.25, 0.25]) + [0.25, -0.25], # Right
    u -> Tor(u * (7π / 4 - 9π / 4) + 9π / 4, 0.0), # Top
    n, n)

joints = ((Joint(2,FaADE.Right),),(Joint(1,FaADE.Left),),)

Dom = GridMultiBlock((D1,D2),joints)



# using GLMakie
using CairoMakie

f = Figure(size=(800,600), fontsize=20)
ax = Axis(f[1,1], xlabel="x", ylabel="y")

g1 = wireframe!(ax,D1.gridx,D1.gridy,zeros(size(D1)))
g2 = wireframe!(ax,D2.gridx,D2.gridy,zeros(size(D2)))


# xlims!(ax, (-2, 1.3))
# ylims!(ax, (-0.4, 1.1))

axleg = axislegend(ax, [g1,g2], ["Block 1","Block 2"], position=:lt)
for i in 1:2
    axleg.entrygroups[][1][2][i].elements[1].attributes[:markersize] = Observable(20)
end
notify(axleg.entrygroups)




save("./figs/Domain_doubleseashell.pdf", f, px_per_unit=8)
