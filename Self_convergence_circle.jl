#=
    Single island error in a circle
=#
using LinearAlgebra
# using CubicHermiteSpline
using Surrogates
using FaADE


using JLD2
using DataFrames
using CSV



# Save the convergence rates?
save_rates          = true

# Run the converging cases?
run_cases           = true
save_cases          = true

# Run the reference solution?
compute_reference   = true
save_reference      = true

# Plot and save?
plot                = true

# Which form?
rejection_form      = false




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


dilation = 0.1 #dilation of the interior domain
g0b = 0.25
α = 0.1

if α == 0.0
    packfn(x) = x
else
    island_location = (0.7 - (g0b + 0.25 * dilation))/(1.0 - (g0b + 0.25 * dilation))
    Dp(x,α,B,B0,B1) = B + α * sinh( asinh((B1-B)/α)*x + asinh((B0-B)/α)*(1-x) )
    packfn(x) = Dp(x, α, island_location, 0.0, 1.0)
end



@show dilation, g0b, α, order
if rejection_form
    dirname = string("Paper2/data/Convergence/SingleIslandSelf_rejection_$(order)_$(EXP)_$(dilation)_$(g0b)_$(α)")
else
    dirname = string("Paper2/data/Convergence/SingleIslandSelf_$(order)_$(EXP)_$(dilation)_$(g0b)_$(α)")
end
@show dirname






### Initial and forcing condition setup

u₀(x,y) = 0.0
F(X,t) = begin
    x,y = X
    4 * (1 - (x^2 + y^2))^8
end
# DIRICHLET
Bxy(X,t) = 0.0   #Boundary condition x=0


### Magnetic field setup

δ = 0.05
rs = 0.7
function B(X,x::Array{Float64},params,t)
    X[1] = δ*x[1]*(1-x[1])*sin(x[2])#/bn
    X[2] = (2x[1] - 2*rs + δ*(1-x[1])*cos(x[2]) - δ*x[1]*cos(x[2]))#/bn
end
dH(X,x,params,t) = B(X,x,params,t)

function MagField(x,t)
    # Get the point in the ψ,θ coordinates
    ψ = sqrt(x[1]^2 + x[2]^2)
    θ = atan(x[2],x[1])

    # Compute the field in these coordinates
    X1 = δ*ψ*(1-ψ)*sin(θ) # B_ψ
    X2 = (2ψ - 2*rₛ + δ*(1-ψ)*cos(θ) - δ*ψ*cos(θ)) # B_θ
    X3 = 0.0

    return [X1,X2,X3]
end

XtoB(x,y) = [sqrt(x^2 + y^2), atan(y,x)]
BtoX(r,θ) = [r*cos(θ), r*sin(θ)]

intercept(u,x,y,t) = begin
    if (sqrt(x^2 + y^2) ≈ 1.0)
        tmp = zero(typeof(x))
        return tmp
    else
        return u
    end
end








function sol(nx,ny,order,packfn,t_f;dt=nothing,MagField=nothing)

    if isnothing(MagField)
        gridoptions = Dict("coords"=>(XtoB,BtoX), "xbound"=>[0.0,1.0], "ybound"=>[0.0,2π])
    else
        gridoptions = Dict("coords"=>(XtoB,BtoX), "xbound"=>[0.0,1.0], "ybound"=>[0.0,2π], "B"=>MagField)
    end
    interpotions = Dict("interpolant"=>:chs,"intercept"=>intercept)

    Dom = FaADE.Grid.squared_circle([(nx,ny) for _ in 1:5], radial_packing=packfn)


    Dr = SAT_Dirichlet(Bxy, Dom.Grids[2].Δy, Up, order, Dom.Grids[2].Δx, :Curvilinear) # Block 2 BCs
    Du = SAT_Dirichlet(Bxy, Dom.Grids[3].Δy, Up, order, Dom.Grids[3].Δx, :Curvilinear) # Block 3 BCs
    Dl = SAT_Dirichlet(Bxy, Dom.Grids[4].Δy, Up, order, Dom.Grids[4].Δx, :Curvilinear) # Block 4 BCs
    Dd = SAT_Dirichlet(Bxy, Dom.Grids[5].Δy, Up, order, Dom.Grids[5].Δx, :Curvilinear) # Block 5 BCs

    BD = Dict(2 => (Dr,), 3 => (Du,), 4 => (Dl,), 5 => (Dd,))


    gdata = construct_grid(dH,Dom,[-2.0π,2.0π],gridoptions=gridoptions)
    PData = FaADE.ParallelOperator.ParallelMultiBlock(gdata,Dom,order,κ=K_para,interpopts=interpotions)


    # Build PDE problem
    P       = Problem2D(order,u₀,K_perp,K_perp,Dom,BD,source=F,parallel=PData)

    # t_f = 1e-2
    # @show Δt = 1.0e-4
    if isnothing(dt)
        @show Δt = 1e-2 * Dom.Grids[1].Δx^2
    else
        Δt = dt
    end
    nf = round(t_f/Δt)
    Δt = t_f/nf

    solve(P,Dom,Δt,1.1Δt)
    soln = solve(P,Dom,Δt,t_f)

    return soln
end





# Interpolate the reference solution back onto the low resolution grid
function reconstruct_soln(Dom,refinterp)
    soln = [zeros(Dom.Grids[I].nx,Dom.Grids[I].ny) for I in 1:5]
    for I in eachindex(Dom.Grids)
        for J in eachindex(Dom.Grids[I].gridx)
            # soln[I][J] = evaluate(refinterp[I],[Dom.Grids[I][J]...])[1]
            soln[I][J] = refinterp[I](tuple(Dom.Grids[I][J]...))
        end
    end
    return soln
end


# Compute the relative error between a low resolution and reconstructed solution
function compute_errors(ublock,refu,Dom)
    tmpnume = 0.0
    tmpdenom = 0.0
    for I in eachindex(ublock) #for each grid

        H = FaADE.Derivatives.innerH(Dom.Grids[I].Δx, Dom.Grids[I].Δy, Dom.Grids[I].nx, Dom.Grids[I].ny, order)

        tmp = ublock[I] .- refu[I]
        tmpnume += H(tmp, Dom.Grids[I].J, tmp)
        tmpdenom += H(refu[I], Dom.Grids[I].J, refu[I]) # (u_ref, H
    end
    return sqrt(tmpnume),sqrt(tmpnume)/sqrt(tmpdenom)
end







println("Start")



if run_cases
    if rejection_form
        solns = [sol(N,N,order,packfn,t_f,MagField=MagField) for N in resolutions]
    else
        solns = [sol(N,N,order,packfn,t_f) for N in resolutions]
    end
    solns_u = [soln.u[2] for soln in solns]

    if save_cases
        try
            mkdir(dirname)
        catch
        end
        jldsave(string(dirname,"/runs");
            solns_u = solns_u
            )
    end
else
    file = jldopen(string(dirname,"/runs"))
    solns_u = file["solns_u"]
end


N = [(n,n) for n in resolutions]
Doms = [FaADE.Grid.squared_circle([n for _ in 1:5], radial_packing=packfn) for n in N]



if compute_reference
    println("reference")
    nx = ny = reference_resolution

    if rejection_form
        reference_soln = sol(nx,ny,4,packfn,t_f,dt=reference_dt,MagField=MagField)
    else
        reference_soln = sol(nx,ny,4,packfn,t_f,dt=reference_dt)
    end
    reference_u = reference_soln.u[2]

    reference_Dom = FaADE.Grid.squared_circle([(nx,ny) for _ in 1:5], radial_packing=packfn)

    if save_reference
        jldsave(string(dirname,"/ref"); reference_u=reference_u)
    end
else
    nx = ny = reference_resolution
    reference_Dom = FaADE.Grid.squared_circle([(nx,ny) for _ in 1:5], radial_packing=packfn)
    @load string(dirname,"/ref") reference_u
end




refinterp = []
for I in eachindex(reference_Dom.Grids)
    grid = tuple.(reference_Dom.Grids[I].gridx[:], reference_Dom.Grids[I].gridy[:])
    push!(refinterp,RadialBasis(grid,reference_u[I],minimum(reference_u[I]),maximum(reference_u[I]),rad=cubicRadial()))
end



reconstructed_solns = [reconstruct_soln(Dom,refinterp) for Dom in Doms]

errors = [compute_errors(soln,reconstructed,Dom) for (soln,reconstructed,Dom) in zip(solns_u, reconstructed_solns,Doms)]

absolute_errors = first.(errors)
relative_errors = last.(errors)


grids = collect(resolutions)*5

println("Completed run for order $(order), with rejection form $(rejection_form).")
@show absolute_errors
@show relative_errors
@show log2.(relative_errors[1:end-1]./relative_errors[2:end]) ./ log2.( grids[2:end] ./ grids[1:end-1] )



if save_rates
    df = DataFrame(N=grids,relative_errors=relative_errors)
    CSV.write("Paper2/data/Convergence/SingleIslandSelfConvergence_O$(order)_EXP$(EXP)_gamma$(dilation)_g0b$(g0b)_alpha$(α).csv",df)
end



if plot
    using GLMakie
    # using LaTeXStrings

    # Solution across the separatrix
    g = Figure()
    axg = Axis(g[1,1], xlabel="x", ylabel="u(x,y)")
    # lines!(axg, soln4.u[2][1][1,:])
    for I in eachindex(solns_u)
        local nx = Doms[I].Grids[4].nx
        lines!(axg,
            Doms[I].Grids[4].gridx[floor(Int,nx/2 + 1),:],
            solns_u[I][4][floor(Int,nx/2 + 1),:],
            label=string(nx))
    end


    h = Figure()
    axh = Axis(h[1,1], xlabel="Resolution", ylabel="Relative error", yscale=log10, xscale=log10)
    scatterlines!(axh, grids, relative_errors)

    h_conv_1 = lines!(axh, grids, grids.^-1, label="O(1)")
    # h_conv_2 = lines!(axh, grids, grids.^-2, label="O(2)")
    # h_conv_4 = lines!(axh, grids, grids.^-4, label="O(4)")

    axislegend(axh, [h_conv_1], ["O(1)"])



    f = Figure(); axf = Axis3(f[1,1])
    for I in 1:5
        surface!(axf, reference_Dom.Grids[I].gridx, reference_Dom.Grids[I].gridy, reference_u[I])
        wireframe!(axf, reference_Dom.Grids[I].gridx, reference_Dom.Grids[I].gridy, reference_u[I])
    end


    save(string(dirname,"/solution.png"),g)
    save(string(dirname,"/convergence.png"),h)
    save(string(dirname,"/reference_solution.png"),f)

end
