#=
    MMS for the circular domain
=#

using LinearAlgebra
using JLD2
using Revise
using FaADE



TestDirichlet = true
SaveTests = true

# Generates the exact MMS solution
function generate_MMS(MMS::Function, grid::GridMultiBlock, t::Float64)
    u_MMS = [zeros(size(grid.Grids[I])) for I in eachindex(grid.Grids)]

    for I in eachindex(grid.Grids)
        for J in eachindex(grid.Grids[I])
            u_MMS[I][J] = MMS(grid.Grids[I][J]..., t)
        end
    end

    return u_MMS
end




function comp_MMS(npts,
    Bxy, BType,
    F, ũ, ũ₀, order;
    dt_scale=1.0, t_f=0.1, kx=1.0, ky=kx, θ=1.0, dilation=0.0, g0b=0.25)

    comp_soln = []
    MMS_soln = []
    grids = []
    # relerr = []
    relerr = zeros(length(npts))

    # Loop
    # for n in npts
    for I in eachindex(npts)
        n = npts[I]

        nx = ny = n

        D1 = Grid2D(
            u -> [-g0b, -g0b] + u * ([g0b, -g0b] - [-g0b, -g0b]) + u * (1 - u) * [0.0, -dilation],
            v -> [-g0b, -g0b] + v * ([-g0b, g0b] - [-g0b, -g0b]) + v * (1 - v) * [-dilation, 0.0],
            v -> [g0b, -g0b] + v * ([g0b, g0b] - [g0b, -g0b]) + v * (1 - v) * [dilation, 0.0],
            u -> [-g0b, g0b] + u * ([g0b, g0b] - [-g0b, g0b]) + u * (1 - u) * [0.0, dilation],
            nx, ny
        )

        T = FaADE.Grid.Torus([1.0], [1.0], [1], [0])

        # Right domain
        D2 = Grid2D(#u->[g0b, -u*0.5 + g0b], # Bottom
            u -> [g0b, g0b] + u * ([g0b, -g0b] - [g0b, g0b]) + u * (1 - u) * [dilation, 0.0],
            v -> v * (T(π / 4, 0.0) - [g0b, g0b]) + [g0b, g0b], # Left
            v -> v * (T(7π / 4, 0.0) + [-g0b, g0b]) + [g0b, -g0b], # Right
            u -> T(u * (7π / 4 - 9π / 4) + 9π / 4, 0.0), # Top
            nx, ny)

        # Top domain
        D3 = Grid2D(#u->[u*0.5 - g0b, g0b], # Bottom
            u -> [-g0b, g0b] + u * ([g0b, g0b] - [-g0b, g0b]) + u * (1 - u) * [0.0, dilation],
            v -> v * (T(3π / 4, 0.0) + [g0b, -g0b]) + [-g0b, g0b], # Left
            v -> v * (T(π / 4, 0.0) - [g0b, g0b]) + [g0b, g0b], # Right
            u -> T(u * (π / 4 - 3π / 4) + 3π / 4, 0.0), # Top
            nx, ny)

        # Left domain
        D4 = Grid2D(#u->[-g0b,u*0.5 - g0b],
            u -> [-g0b, -g0b] + u * ([-g0b, g0b] - [-g0b, -g0b]) + u * (1 - u) * [-dilation, 0.0],
            v -> v * (T(5π / 4, 0.0) - [-g0b, -g0b]) + [-g0b, -g0b],
            v -> v * (T(3π / 4, 0.0) - [-g0b, g0b]) + [-g0b, g0b],
            u -> T(u * (3π / 4 - 5π / 4) + 5π / 4, 0.0),
            nx, ny)

        # Bottom domain
        D5 = Grid2D(#u->[-u*0.5 + g0b, -g0b],
            u -> [-g0b, -g0b] + u * ([g0b, -g0b] - [-g0b, -g0b]) + u * (1 - u) * [0.0, -dilation],
            v -> v * (T(7π / 4, 0.0) - [g0b, -g0b]) + [g0b, -g0b],
            v -> v * (T(5π / 4, 0.0) - [-g0b, -g0b]) + [-g0b, -g0b],
            u -> T(u * (5π / 4 - 7π / 4) + 7π / 4, 0.0),
            nx, ny)


        joints = ((Joint(2, Right), Joint(3, Up), Joint(4, Left), Joint(5, Down)),
            (Joint(1, Down), Joint(3, Left), Joint(5, Right)),
            (Joint(1, Down), Joint(4, Left), Joint(2, Right)),
            (Joint(1, Down), Joint(5, Left), Joint(3, Right)),
            (Joint(1, Down), Joint(2, Left), Joint(4, Right)))


        Dom = GridMultiBlock((D1, D2, D3, D4, D5), joints)

        if BType == Dirichlet


            Dr = FaADE.SATs.SAT_Dirichlet(Bxy, D2.Δy, Up, order, D2.Δx, :Curvilinear) # Block 2 BCs
            Du = FaADE.SATs.SAT_Dirichlet(Bxy, D3.Δy, Up, order, D3.Δx, :Curvilinear) # Block 3 BCs
            Dl = FaADE.SATs.SAT_Dirichlet(Bxy, D4.Δy, Up, order, D4.Δx, :Curvilinear) # Block 4 BCs
            Dd = FaADE.SATs.SAT_Dirichlet(Bxy, D5.Δy, Up, order, D5.Δx, :Curvilinear) # Block 5 BCs

            BD = Dict(2 => (Dr,), 3 => (Du,), 4 => (Dl,), 5 => (Dd,))


        elseif BType == Neumann
        end


        Δt = dt_scale * min(D1.Δx, D1.Δy)^2
        t_f = 1e-1
        nt = round(t_f / Δt)
        Δt = t_f / nt

        P = Problem2D(order, ũ₀, K, K, Dom, BD, source=F)

        println("Solving n=", n, " case with Δt=", Δt)
        soln = solve(P, Dom, Δt, t_f)

        u_MMS = generate_MMS(ũ, Dom, soln.t[2])

        push!(comp_soln, soln)
        push!(grids, Dom)
        push!(MMS_soln, u_MMS)

        tmpnume = 0.0
        tmpdeno = 0.0

        for J in eachindex(Dom.Grids)
            H = FaADE.Derivatives.innerH(Dom.Grids[J].Δx, Dom.Grids[J].Δy, Dom.Grids[J].nx, Dom.Grids[J].ny, order)
            v = soln.u[2][J] .- u_MMS[J]
            tmpnume += H(v, Dom.Grids[J].J, v)
            tmpdeno += H(u_MMS[J], Dom.Grids[J].J, u_MMS[J])
        end
        relerr[I] = sqrt(tmpnume) / sqrt(tmpdeno)

    end

    @show conv_rate = log.(relerr[1:end-1] ./ relerr[2:end]) ./ log.((1 ./ (npts[1:end-1] .- 1)) ./ (1 ./ (npts[2:end] .- 1)))
    # conv_rate = log2.(relerr[1][1:end-1] ./ relerr[1][2:end], relerr[2][1:end-1] ./ relerr[2][2:end])
    # @show conv_rate = log2.(relerr[1:end-1] ./ relerr[2:end])

    return (comp_soln=comp_soln, MMS_soln=MMS_soln, grids=grids, relerr=relerr, conv_rate=conv_rate, npts=npts)
end



########################################################################################
########################################################################################
# CALLING
########################################################################################
########################################################################################

order = [2, 4]
npts = [11, 21, 31, 41, 51, 61]
dilation = [0.0, 0.1]
inner_bound = [0.25]

θ = 0.5

# Solution
ũ(x, y, t;
    ωt=1.0,
    ωx=1.0, cx=0.0,
    ωy=1.0, cy=0.0) = cos(2π * ωt * t) * sin(2π * x * ωx + cx) * sin(2π * y * ωy + cy)

# Initial condition
ũ₀(x, y;
    ωx=1.0, cx=0.0,
    ωy=1.0, cy=0.0) = sin(2π * ωx * x + cx) * sin(2π * ωy * y + cy)


K = 1.0

F(x, y, t;
    ωt=1.0,
    ωx=1.0, cx=0.0,
    ωy=1.0, cy=0.0,
    K=1.0) =
    -2π * ωt * sin(2π * ωt * t) * sin(2π * x * ωx + cx) * sin(2π * y * ωy + cy) +
    K * 4π^2 * ωx^2 * cos(2π * ωt * t) * sin(2π * x * ωx + cx) * sin(2π * y * ωy + cy) +
    K * 4π^2 * ωy^2 * cos(2π * ωt * t) * sin(2π * x * ωx + cx) * sin(2π * y * ωy + cy)


if TestDirichlet
    println("=====")
    println("Dirichlet")
    cx = 0.0
    cy = 0.0
    ωx = 5.5
    ωy = 6.0
    ωt = 3.0

    println("ωx=", ωx, "  ωy=", ωy, ",  cx=", cx, ",  cy=", cy, ", ωt=", ωt, " θ=", θ)

    analytic(x, y, t) = ũ(x, y, t, ωt=ωt, ωx=ωx, cx=cx, ωy=ωy, cy=cy)
    IC(x, y) = ũ₀(x, y, ωx=ωx, cx=cx, ωy=ωy, cy=cy)
    FD(X, t) = F(X[1], X[2], t, ωt=ωt, ωx=ωx, cx=cx, ωy=ωy, cy=cy, K=K)

    Bũ(X, t) = cos(2π * ωt * t) * sin(2π * ωx * X[1] + cx) * sin(2π * X[2] * ωy + cy) #Boundary condition x=0

    Dirichlet_MMS = Dict()
    for ord in order
        Dirichlet_MMS[ord] = Dict()

        for dil in dilation
            Dirichlet_MMS[ord][dil] = Dict()

            for g0b in inner_bound
                println("Order=", ord, " Dilation=", dil, " boundary=", g0b)
                Dirichlet_MMS[ord][dil][g0b] = comp_MMS(npts,
                    Bũ, Dirichlet,
                    FD, analytic, IC, ord,
                    kx=K, ky=K, θ=θ, dilation=dil, g0b=g0b)
            end
        end
    end

    println("=====")
end





if SaveTests

    try
        mkdir("Paper2/MMSData")
    catch
    end

    for ord in order, dil in dilation, g0b in inner_bound
        tmp = [Dirichlet_MMS[ord][dil][g0b].comp_soln[I].u for I in eachindex(npts)]
        tmp_mms = [Dirichlet_MMS[ord][dil][g0b].MMS_soln[I] for I in eachindex(npts)]
        jldsave(string("Paper2/MMSData/MMS_$(ord)_$(g0b)_$(dil)");
            solns=tmp,
            mms_solns=tmp_mms,
            npts=npts
        )
    end


    nameappend = string("conv")

    using DataFrames, CSV
    df = DataFrame()


    headervals = [(ord, dil, g0b) for ord in order, dil in dilation, g0b in inner_bound][:]
    header = vcat(["N "], ["$(ord)_$(dil)_$(g0b)" for ord in order, dil in dilation, g0b in inner_bound][:]...)

    relerrdict = Dict()
    relerrdict["N"] = npts

    for val in headervals
        relerrdict[replace(string(val)[2:end-1], "," => "_", " " => "")] = Dirichlet_MMS[val[1]][val[2]][val[3]].relerr
    end

    df = DataFrame(relerrdict)

    select!(df, :N, Not(:N)) #Ensure N is the first column


    CSV.write("./Paper2/data/MMS_Circle.csv", df)

end
