#=
    Pictures of the circle domain with a square domain and a bowed domain.
=#
# using GLMakie
using CairoMakie
using FaADE




function buildgrid(nx, ny, dialation)

    D1 = Grid2D(
        u -> [-0.25, -0.25] + u * ([0.25, -0.25] - [-0.25, -0.25]) + u * (1 - u) * [0.0, -dialation],
        v -> [-0.25, -0.25] + v * ([-0.25, 0.25] - [-0.25, -0.25]) + v * (1 - v) * [-dialation, 0.0],
        v -> [0.25, -0.25] + v * ([0.25, 0.25] - [0.25, -0.25]) + v * (1 - v) * [dialation, 0.0],
        u -> [-0.25, 0.25] + u * ([0.25, 0.25] - [-0.25, 0.25]) + u * (1 - u) * [0.0, dialation],
        nx, ny
    )

    Tor = FaADE.Grid.Torus([1.0], [1.0], [1], [0])

    # Right domain
    D2 = Grid2D(#u->[0.25, -u*0.5 + 0.25], # Bottom
        u -> [0.25, 0.25] + u * ([0.25, -0.25] - [0.25, 0.25]) + u * (1 - u) * [dialation, 0.0],
        v -> v * (Tor(π / 4, 0.0) - [0.25, 0.25]) + [0.25, 0.25], # Left
        v -> v * (Tor(7π / 4, 0.0) + [-0.25, 0.25]) + [0.25, -0.25], # Right
        u -> Tor(u * (7π / 4 - 9π / 4) + 9π / 4, 0.0), # Top
        nx, ny)

    # Top domain
    D3 = Grid2D(#u->[u*0.5 - 0.25, 0.25], # Bottom
        u -> [-0.25, 0.25] + u * ([0.25, 0.25] - [-0.25, 0.25]) + u * (1 - u) * [0.0, dialation],
        v -> v * (Tor(3π / 4, 0.0) + [0.25, -0.25]) + [-0.25, 0.25], # Left
        v -> v * (Tor(π / 4, 0.0) - [0.25, 0.25]) + [0.25, 0.25], # Right
        u -> Tor(u * (π / 4 - 3π / 4) + 3π / 4, 0.0), # Top
        nx, ny)

    # Left domain
    D4 = Grid2D(#u->[-0.25,u*0.5 - 0.25],
        u -> [-0.25, -0.25] + u * ([-0.25, 0.25] - [-0.25, -0.25]) + u * (1 - u) * [-dialation, 0.0],
        v -> v * (Tor(5π / 4, 0.0) - [-0.25, -0.25]) + [-0.25, -0.25],
        v -> v * (Tor(3π / 4, 0.0) - [-0.25, 0.25]) + [-0.25, 0.25],
        u -> Tor(u * (3π / 4 - 5π / 4) + 5π / 4, 0.0),
        nx, ny)

    # Bottom domain
    D5 = Grid2D(#u->[-u*0.5 + 0.25, -0.25],
        u -> [-0.25, -0.25] + u * ([0.25, -0.25] - [-0.25, -0.25]) + u * (1 - u) * [0.0, -dialation],
        v -> v * (Tor(7π / 4, 0.0) - [0.25, -0.25]) + [0.25, -0.25],
        v -> v * (Tor(5π / 4, 0.0) - [-0.25, -0.25]) + [-0.25, -0.25],
        u -> Tor(u * (5π / 4 - 7π / 4) + 7π / 4, 0.0),
        nx, ny)


    joints = ((Joint(2, FaADE.Right), Joint(3, FaADE.Up), Joint(4, FaADE.Left), Joint(5, FaADE.Down)),
        (Joint(1, FaADE.Down), Joint(3, FaADE.Left), Joint(5, FaADE.Right)),
        (Joint(1, FaADE.Down), Joint(4, FaADE.Left), Joint(2, FaADE.Right)),
        (Joint(1, FaADE.Down), Joint(5, FaADE.Left), Joint(3, FaADE.Right)),
        (Joint(1, FaADE.Down), Joint(2, FaADE.Left), Joint(4, FaADE.Right)))


    Dom = GridMultiBlock((D1, D2, D3, D4, D5), joints)

    return Dom
end









#=== Generate figure for the circle ===#

nx = ny = 11


Dom0 = buildgrid(nx, ny, 0.0)
Dom2 = buildgrid(nx, ny, 0.1)


gridfig = Figure(size=(1200, 600), fontsize=25)

gridfig_gg = gridfig[1, 1] = GridLayout()

gridfig_ax1 = Axis(gridfig_gg[1, 1], ylabel="Z", xlabel="R")
gridfig_ax2 = Axis(gridfig_gg[1, 2], xlabel="R")

dg = []

# Outline of exterior volumes
for I in 1:5
    wireframe!(gridfig_ax1, Dom0.Grids[I].gridx, Dom0.Grids[I].gridy, zeros(nx, ny))
    push!(dg, wireframe!(gridfig_ax2, Dom2.Grids[I].gridx, Dom2.Grids[I].gridy, zeros(nx, ny)))
end

# Outline of exterior volumes
for I in 2:5
    lines!(gridfig_ax1, Dom0.Grids[I].gridx[:, 1], Dom0.Grids[I].gridy[:, 1], color=(:black, 0.7))
    lines!(gridfig_ax1, Dom0.Grids[I].gridx[end, :], Dom0.Grids[I].gridy[end, :], color=(:black, 0.7))
    lines!(gridfig_ax1, Dom0.Grids[I].gridx[:, end], Dom0.Grids[I].gridy[:, end], color=(:black, 0.7))
end

# Outline of exterior volumes
for I in 2:5
    lines!(gridfig_ax2, Dom2.Grids[I].gridx[:, 1], Dom2.Grids[I].gridy[:, 1], color=(:black, 0.7))
    lines!(gridfig_ax2, Dom2.Grids[I].gridx[end, :], Dom2.Grids[I].gridy[end, :], color=(:black, 0.7))
    lines!(gridfig_ax2, Dom2.Grids[I].gridx[:, end], Dom2.Grids[I].gridy[:, end], color=(:black, 0.7))
end


linkyaxes!(gridfig_ax1, gridfig_ax2)

hideydecorations!(gridfig_ax2, grid=false)


Legend(gridfig[1, 2], dg[:], [L"D_1", L"D_2", L"D_3", L"D_4", L"D_5"])

text!(gridfig_ax1, 1, 1, align=(:right, :top), text=L"\gamma=0")
text!(gridfig_ax2, 1, 1, align=(:right, :top), text=L"\gamma=0.1")



# Safe figure to ./data/ folder
save("./data/F_Domain_bowed.pdf", gridfig)
