using GLMakie

using SPECreader


fname = "Paper2/data/solovev_pert_a_pap.sp.h5"


speceq = SPECEquilibrium(fname)

# specpoin = ReadPoincare(fname)


poincare = SPECreader.construct_poincare(speceq, N_trajs=200, N_orbs=500)

boundary = [SPECreader.get_RZ(1.0, t, 0.0, speceq, 2) for t in range(0.0, 2π, 101)]

f = Figure(size=(5, 3) .* 150)
axf = Axis(f[1, 1])
scatter!(axf, poincare.ψ, poincare.θ, markersize=2.0)
lines!(axf, first.(boundary), last.(boundary), color=:black)

xlims!(minimum(first.(boundary)) - 1e-2, maximum(first.(boundary)) + 1e-2)
# ylims!(minimum(last.(boundary)) - 1e-2, maximum(last.(boundary)) + 1e-2)
ylims!(0.0, maximum(last.(boundary)) + 1e-2)


hidedecorations!(axf)


f

using CairoMakie
CairoMakie.activate!()
save("./Paper2/data/IntroPlot.png", f, px_per_unit=4)
