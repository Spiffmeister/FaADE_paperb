#=
    Pictures of the circle domain with a square domain and a bowed domain.
=#

using LinearAlgebra
using Revise
using FaADE



#====== New solver 4 volume ======#
println("Curvilinear volume")



function buildgrid(nx,ny,dialation)

    D1 = Grid2D(
        u -> [-0.25,-0.25] + u*([0.25,-0.25] - [-0.25,-0.25]) + u*(1-u)*[0.0,-dialation],
        v -> [-0.25,-0.25] + v*([-0.25,0.25] - [-0.25,-0.25]) + v*(1-v)*[-dialation,0.0],
        v -> [0.25,-0.25] + v*([0.25,0.25] - [0.25,-0.25]) + v*(1-v)*[dialation,0.0],
        u -> [-0.25,0.25] + u*([0.25,0.25] - [-0.25,0.25]) + u*(1-u)*[0.0,dialation],
        nx,ny
    )

    Tor = FaADE.Grid.Torus([1.0],[1.0],[1],[0])

    # Right domain
    D2 = Grid2D(#u->[0.25, -u*0.5 + 0.25], # Bottom
                u->[0.25,0.25] + u*([0.25,-0.25] - [0.25,0.25]) + u*(1-u)*[dialation,0.0],
                v->v*(Tor(π/4,0.0) - [0.25,0.25]) + [0.25,0.25], # Left
                v->v*(Tor(7π/4,0.0) + [-0.25, 0.25]) + [0.25, -0.25], # Right
                u->Tor(u*(7π/4 - 9π/4) + 9π/4,0.0), # Top
                nx,ny)

    # Top domain
    D3 = Grid2D(#u->[u*0.5 - 0.25, 0.25], # Bottom
                u->[-0.25,0.25] + u*([0.25,0.25] - [-0.25,0.25]) + u*(1-u)*[0.0,dialation],
                v->v*(Tor(3π/4,0.0) + [0.25,-0.25]) + [-0.25,0.25], # Left
                v->v*(Tor(π/4,0.0) - [0.25,0.25]) + [0.25,0.25], # Right
                u->Tor(u*(π/4 - 3π/4) + 3π/4,0.0), # Top
                nx,ny)

    # Left domain
    D4 = Grid2D(#u->[-0.25,u*0.5 - 0.25],
                u->[-0.25,-0.25] + u*([-0.25,0.25] - [-0.25,-0.25]) + u*(1-u)*[-dialation,0.0],
                v->v*(Tor(5π/4,0.0) - [-0.25, -0.25]) + [-0.25, -0.25],
                v->v*(Tor(3π/4,0.0) - [-0.25,0.25]) + [-0.25,0.25],
                u->Tor(u*(3π/4 - 5π/4) + 5π/4,0.0),
                nx,ny)

    # Bottom domain
    D5 = Grid2D(#u->[-u*0.5 + 0.25, -0.25],
                u->[-0.25,-0.25] + u*([0.25,-0.25] - [-0.25,-0.25]) + u*(1-u)*[0.0,-dialation],
                v->v*(Tor(7π/4,0.0) - [0.25,-0.25]) + [0.25, -0.25],
                v->v*(Tor(5π/4,0.0) - [-0.25,-0.25]) + [-0.25, -0.25],
                u->Tor(u*(5π/4 - 7π/4) + 7π/4, 0.0),
                nx,ny)


    joints = ((Joint(2,FaADE.Right),Joint(3,FaADE.Up),Joint(4,FaADE.Left),Joint(5,FaADE.Down)),
                (Joint(1,FaADE.Down),Joint(3,FaADE.Left),Joint(5,FaADE.Right)),
                (Joint(1,FaADE.Down),Joint(4,FaADE.Left),Joint(2,FaADE.Right)),
                (Joint(1,FaADE.Down),Joint(5,FaADE.Left),Joint(3,FaADE.Right)),
                (Joint(1,FaADE.Down),Joint(2,FaADE.Left),Joint(4,FaADE.Right)))


    Dom = GridMultiBlock((D1,D2,D3,D4,D5),joints)

    return Dom
end

#=
δ = 0.05
rs = 0.5
function B(X,x::Array{Float64},params,t)
    # X[1] = -x[1] * δ * (-x[1]^4 + 1) * sin(x[2])
    # X[2] = -2x[1] + 2rs - 2δ * x[1] * (-x[1]^4 + 1) * cos(x[2]) + 4δ * x[1]^5 * cos(x[2])

    X[1] = δ*x[1]*(1-x[1])*sin(x[2])#/bn
    X[2] = (2x[1] - 2*rs + δ*(1-x[1])*cos(x[2]) - δ*x[1]*cos(x[2]))#/bn
    # time dependent
    # X[1] = -x[1] * δ * (-r^4 + 1) * sin(x[2]) * t
    # X[1] = (-2x[1] + 2rs - 2δ * x[1] * (-x[1]^4 + 1) * cos(x[2]) + 4δ * r^5 * cos(x[2])) * t
end
dH(X,x,params,t) = B(X,x,params,t)

XtoB(x,y) = [sqrt(x^2 + y^2), atan(y,x)]
BtoX(r,θ) = [r*cos(θ), r*sin(θ)]
gridoptions = Dict("coords"=>(XtoB,BtoX), "xbound"=>[0.0,1.0], "ybound"=>[0.0,2π], "remapping"=>:bilinear)

gdata = construct_grid(dH,Dom,[-2.0π,2.0π],gridoptions=gridoptions)

# gdata = remap_grid(gdata,interpmode=:idw)

PData = FaADE.ParallelOperator.ParallelMultiBlock(gdata,Dom,order,κ=1.0e6)
=#







# using GLMakie
using CairoMakie



nx = ny = 11


Dom0 = buildgrid(nx,ny,0.0)
Dom2 = buildgrid(nx,ny,0.1)
# Dom4 = buildgrid(nx,ny,0.4)


gridfig = Figure(size=(1200,600), fontsize=20)

gridfig_gg = gridfig[1,1] = GridLayout()

gridfig_ax1 = Axis(gridfig_gg[1,1], ylabel="Z", xlabel="R")
gridfig_ax2 = Axis(gridfig_gg[1,2],             xlabel="R")
# gridfig_ax3 = Axis(gridfig_gg[1,3],             xlabel="R")

dg = []

for I in 1:5
    wireframe!(gridfig_ax1,Dom0.Grids[I].gridx,Dom0.Grids[I].gridy,zeros(nx,ny))
    push!(dg, wireframe!(gridfig_ax2,Dom2.Grids[I].gridx,Dom2.Grids[I].gridy,zeros(nx,ny)))
    # push!(dg,wireframe!(gridfig_ax3,Dom4.Grids[I].gridx,Dom4.Grids[I].gridy,zeros(nx,ny)))
end

for I in 2:5
    # lines!(gridfig_ax1,DomSquare.Grids[I].gridx[1,:],   DomSquare.Grids[I].gridy[1,:],  color=(:black,0.5))
    lines!(gridfig_ax1,Dom0.Grids[I].gridx[:,1],   Dom0.Grids[I].gridy[:,1],  color=(:black,0.7))
    lines!(gridfig_ax1,Dom0.Grids[I].gridx[end,:], Dom0.Grids[I].gridy[end,:],color=(:black,0.7))
    lines!(gridfig_ax1,Dom0.Grids[I].gridx[:,end], Dom0.Grids[I].gridy[:,end],color=(:black,0.7))
end


for I in 2:5
    # lines!(gridfig_ax2,DomBowed.Grids[I].gridx[1,:],   DomBowed.Grids[I].gridy[1,:],  color=:black)
    lines!(gridfig_ax2,Dom2.Grids[I].gridx[:,1],   Dom2.Grids[I].gridy[:,1],  color=(:black,0.7))
    lines!(gridfig_ax2,Dom2.Grids[I].gridx[end,:], Dom2.Grids[I].gridy[end,:],color=(:black,0.7))
    lines!(gridfig_ax2,Dom2.Grids[I].gridx[:,end], Dom2.Grids[I].gridy[:,end],color=(:black,0.7))
end

# for I in 2:5
    # lines!(gridfig_ax2,DomBowed.Grids[I].gridx[1,:],   DomBowed.Grids[I].gridy[1,:],  color=:black)
    # lines!(gridfig_ax3,Dom4.Grids[I].gridx[:,1],   Dom4.Grids[I].gridy[:,1],  color=(:black,0.7))
    # lines!(gridfig_ax3,Dom4.Grids[I].gridx[end,:], Dom4.Grids[I].gridy[end,:],color=(:black,0.7))
    # lines!(gridfig_ax3,Dom4.Grids[I].gridx[:,end], Dom4.Grids[I].gridy[:,end],color=(:black,0.7))
# end


linkyaxes!(gridfig_ax1,gridfig_ax2)
# linkyaxes!(gridfig_ax1,gridfig_ax2,gridfig_ax3)

hideydecorations!(gridfig_ax2, grid=false)
# hideydecorations!(gridfig_ax3, grid=false)


# axislegend(gridfig_ax3, dg[:], [L"D_1",L"D_2",L"D_3",L"D_4",L"D_5"], position=:lt)
# axislegend(gridfig_ax2, dg[:], [L"D_1",L"D_2",L"D_3",L"D_4",L"D_5"], position=:lt)
Legend(gridfig[1,2], dg[:], [L"D_1",L"D_2",L"D_3",L"D_4",L"D_5"])

# Label(gridfig_gg[1,1,Makie.TopRight()], L"\gamma=0")
text!(gridfig_ax1,1,1,align=(:right,:top),text=L"\gamma=0")
text!(gridfig_ax2,1,1,align=(:right,:top),text=L"\gamma=0.1")
# text!(gridfig_ax3,1,1,align=(:right,:top),text=L"\gamma=0.4")



# gridfig
save("Paper2/data/Domain_bowed.pdf",gridfig)






# jacfig = Figure(size=(1200,500), fontsize=20)
# jacfig_gg = jacfig[1,1] = GridLayout()

# jacfig_ax1 = Axis3(jacfig_gg[1,1], ylabel="Z", xlabel="R")
# jacfig_ax2 = Axis3(jacfig_gg[1,2], ylabel="Z", xlabel="R")

# surface!(jacfig_ax1,DomSquare.Grids[1].gridx,DomSquare.Grids[1].gridy,DomSquare.Grids[1].J)
# surface!(jacfig_ax1,DomSquare.Grids[2].gridx,DomSquare.Grids[2].gridy,DomSquare.Grids[2].J)



# surface!(jacfig_ax2,DomBowed.Grids[1].gridx,DomBowed.Grids[1].gridy,DomBowed.Grids[1].J)
# surface!(jacfig_ax2,DomBowed.Grids[2].gridx,DomBowed.Grids[2].gridy,DomBowed.Grids[2].J)


