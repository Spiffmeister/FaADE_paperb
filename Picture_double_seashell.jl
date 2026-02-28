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

θ₀ = π / 6
θ₁ = 2π - θ₀
θ₂ = 2π + θ₀
D2 = Grid2D(#u->[0.25, -u*0.5 + 0.25], # Bottom
    u -> [0.25, 0.25] + u * ([0.25, -0.25] - [0.25, 0.25]) + u * (1 - u) * [dilation, 0.0],
    v -> v * (Tor(θ₀, 0.0) - [0.25, 0.25]) + [0.25, 0.25], # Left
    v -> v * (Tor(θ₁, 0.0) + [-0.25, 0.25]) + [0.25, -0.25], # Right
    u -> Tor(u * (θ₁ - θ₂) + θ₂, 0.0), # Top
    n, n)

joints = ((Joint(2,FaADE.Right),),(Joint(1,FaADE.Left),),)

Dom = GridMultiBlock((D1,D2),joints)



D21 = Grid2D([0.0, 1.0], [0.0, 1.0], 11, 11)
D22 = Grid2D([0.0, 1.0], [0.0, 1.0], 11, 11)


#### FIGURE


# using GLMakie
using CairoMakie


alignmode_top = 0
alignmode_bottom = 0



f = Figure(size=(1200, 450), fontsize=20, figure_padding=1)
gax1 = f[1, 1] = GridLayout(alignmode=Mixed(top=alignmode_top, bottom=alignmode_bottom))

ax1 = Axis(gax1[1, 1], xlabel="x", ylabel="y", aspect=DataAspect())

g1 = wireframe!(ax1, D1.gridx, D1.gridy, zeros(size(D1)))
g2 = wireframe!(ax1, D2.gridx, D2.gridy, zeros(size(D2)))


# xlims!(ax, (-2, 1.3))
# ylims!(ax, (-0.4, 1.1))

axleg = axislegend(ax1, [g1, g2], ["Block 1", "Block 2"], position=:lt)
for i in 1:2
    axleg.entrygroups[][1][2][i].elements[1].attributes[:markersize] = Observable(20)
end
notify(axleg.entrygroups)






# Place an arrow between the plots
garrow = f[1, 2] = GridLayout(alignmode=Mixed(top=alignmode_top, bottom=50))
garrowax = Axis(garrow[1, 1], aspect=DataAspect(), xautolimitmargin=(0, 0), yautolimitmargin=(0, 0))
poly!(garrowax, Point2f[(0, 0), (0.5, 0), # tail to head
        (0.5, -0.1), (0.8, 0.05), (0.5, 0.2), #head
        (0.5, 0.1), (0, 0.1)], #head to tail
    color=:black)
colsize!(f.layout, 2, Relative(1 / 6))
hidedecorations!(garrowax)
hidespines!(garrowax)




gax2 = f[1, 3] = GridLayout(alignmode=Mixed(top=alignmode_top, bottom=35))

xticks_values = range(0.0, 1.0, 3) |> collect
xtick_strings = fill("", 3)



ax21 = Axis(gax2[1, 1], xlabel="q", title="Block 1", aspect=DataAspect(), xticks=(xticks_values, xtick_strings), xlabelpadding=15.0)
ax22 = Axis(gax2[1, 2], xlabel="q", yaxisposition=:right, ylabel="r", aspect=DataAspect(), title="Block 2", xticks=(xticks_values, xtick_strings), xlabelpadding=15.0)
linkyaxes!(ax21, ax22)
# linkyaxes!(ax21, ax211)

# linkxaxes!(ax21, ax22)

g21 = wireframe!(ax21, D21.gridx, D21.gridy, zeros(size(D21)))
g22 = wireframe!(ax22, D22.gridx, D22.gridy, zeros(size(D22)), color=g2.color)

hideydecorations!(ax21)
colgap!(gax2, 0.0)

xlims!(ax21, low=nothing, high=1.0)
xlims!(ax22, low=0.0, high=nothing)

# Add second axis objects, hide interiors, ensure correct scaling but adding the wireframe but making alpha=0 and then disable all but x-axis labels
ax211 = Axis(gax2[1, 1], aspect=DataAspect(), xticks=([0.1, 0.5, 0.9], string.(xticks_values)), xticksize=0.0)
wireframe!(ax211, D21.gridx, D21.gridy, zeros(size(D21)), alpha=0.0)
xlims!(ax211, low=nothing, high=1.0)
hideydecorations!(ax211)


ax221 = Axis(gax2[1, 2], aspect=DataAspect(), xticks=([0.1, 0.5, 0.9], string.(xticks_values)), xticksize=0.0)
wireframe!(ax221, D22.gridx, D22.gridy, zeros(size(D22)), alpha=0.0)
xlims!(ax221, low=0.0, high=nothing)
hideydecorations!(ax221)


colgap!(f.layout, 1, 5.0)
colgap!(f.layout, 2, 5.0)

resize_to_layout!(f)
# f

save("./Paper2/data/Domain_doubleseashell.pdf", f, px_per_unit=8)
