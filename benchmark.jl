include("ende.jl")
include("endeplots.jl")

using FileIO
using Dates
using UnPack

N = 250

Bmag = 5e-5 # Teslas (kg/As^2)
B_hat = [0.; 0.; 1.]
J = [41.87 0.0 0.0; 0.0 41.87 0.0; 0.0 0.0 6.67] * 1e-3 # kg*m^2
ωmax = 5*π/180
default = Dict(:B_I => B_hat * Bmag,
               :m_max => 0.66,
               :J => J,
               :ωmax => ωmax,
               :n => 1e6,
               :ϵ => 0.0,
               :samplestate => ()->ComponentArray(q = toscalarlast(randquat()), ω = randinunit(3) * ωmax),
               :samplecontrol => ()->[rand()*0.005; rand()*0.03],
               :sampletime => ()->rand()*25,
               :distfn => (x1, x2)->dist(x1, x2, 30),
               :controlfnfact => (x...)->planarpdcontrolfactory(x...),
               :pbias => 0.66,
               :kbias => 25,
               :biasfactor => Inf,
               :pbend => 0.03,
               :kbend => 50,
               :bendfactor => Inf,
               :tmod => 0.1,
               :kmod => nothing,
               :globalbranchimprovement => true,
               :timeout => 300
)

cases = Dict{String,Dict{Symbol,Any}}()

cases["bonsai-pd"] = deepcopy(default)
cases["bonsai-pd"][:controlfnfact] = (x...)->planarpdcontrolfactory(x...)
cases["bonsai-pd"][:samplecontrol] = ()->[rand()*0.005; rand()*0.03]
cases["bonsai-pd"][:timeout] = Inf
cases["bonsai-pd"][:n] = 10000
cases["bonsai-pd"][:ϵ] = 0.01

outputs = Dict{String,Union{Symbol,Function}}()
outputs["time"] = ()->nothing
outputs["cost"] = (path, roadmap, states, controls, prob, distfn)->timecost(controls, path)
outputs["qerror"] = (path, roadmap, states, controls, prob, distfn)->dist(states[path[end]], stategoal(prob), 0)
outputs["werror"] = (path, roadmap, states, controls, prob, distfn)->dist(states[path[end]], stategoal(prob), 1e10)/1e10
outputs["error"] = (path, roadmap, states, controls, prob, distfn)->dist(states[path[end]], stategoal(prob), 30)
outputs["treesize"] = (path, roadmap, states, controls, prob, distfn)->length(states)

out = Dict{String,Dict{String,Vector}}(case => Dict(output => zeros(N) for output in keys(outputs)) for case in keys(cases))

for case in collect(keys(cases))
    vals = cases[case]
    @info "Benchmarking" case vals
    @unpack B_I, m_max, J, ωmax, n, ϵ, samplestate, samplecontrol, sampletime, distfn, controlfnfact, pbias, kbias, biasfactor, pbend, kbend, bendfactor, tmod, kmod, globalbranchimprovement, timeout = vals

    Threads.@threads for i in 1:N
        # prob = AttitudeProblem(Bmag, m_max, J, ωmax, randquat(), Keepout[])
        prob = AttitudeProblem(B_I, m_max, J, ωmax, Quaternion(0.,0.,0.,1.), [Keepout([0.,0.,-1.],[-1.,0.,0.],50*π/180),Keepout(normalize([1.,1.,0.]),[0.,-1.,0.],20*π/180)])
        t0 = time()
        path, roadmap, states, controls = bonsairrt(prob, controlfnfact, samplecontrol, samplestate, sampletime, distfn, n, ϵ, pbias, kbias, biasfactor, pbend, kbend, bendfactor; tmod=tmod, kmod=kmod, globalbranchimprovement=globalbranchimprovement, timeout=timeout)
        telapsed = time()-t0

        for (output, expr) in outputs
            if output == "time"
                out[case][output][i] = telapsed
            else
                out[case][output][i] = expr(path, roadmap, states, controls, prob, distfn)
            end
        end
    end
end
FileIO.save(string("benchmark_",now(),".jld2"),"run",[out,cases])


# Boxplots
groupings = [["rrt-const","rrt-pd"],["bonsai-const","bonsai-pd"]]
# groupings = [["bonsai-pd","bonsaibiasfac-pd","bonsaibendfac-pd","bonsaibiasbendfac-pd"]]
results = ["error","treesize","cost"]
labels = [L"Minimum $d$", "# Nodes", "Maneuver Length [s]"]

colorlist = Makie.wong_colors()
fig = Figure(resolution=(6.5,2).*72.0.*2.5, font = "CMU Serif", fontsize = 18)
for (axi, result) in enumerate(results)
    ax = Axis(fig[1,axi])
    xs = vcat([repeat([g],length(out[trial][result])) for (g,group) in enumerate(groupings) for (t,trial) in enumerate(group) ]...)
    ys = vcat([out[trial][result] for (g,group) in enumerate(groupings) for (t,trial) in enumerate(group) ]...)
    dodge = vcat([repeat([t],length(out[trial][result])) for (g,group) in enumerate(groupings) for (t,trial) in enumerate(group) ]...)
    boxplot!(ax, xs, ys, dodge = dodge, color = [colorlist[d] for d in dodge]) # show_notch = true
    ax.ylabel = labels[axi]
    ax.xticks = ([1,2], ["RRT","Bonsai"])
    if result == "cost"
        ylims!(ax, (0.,1000.))
    end
end

elem_1 = [PolyElement(color = colorlist[1], strokecolor = :black, strokewidth = 1)]
elem_2 = [PolyElement(color = colorlist[2], strokecolor = :black, strokewidth = 1)]
Legend(fig[1, 4], [elem_1, elem_2], ["Const.", "PD"], "Control", titlefont="CMU Serif Bold")

save("output/benchmark_compare.pdf", fig, pt_per_unit=1/2.5)



# Case Compare
fig = Figure(resolution=(6.5,3).*72.0.*2.5, font = "CMU Serif", fontsize = 18)
# fig = Figure(resolution=(1280,1280), fontsize=12)
ax1 = Axis(fig[1,1])
ax2 = Axis(fig[1,2])
ax3 = Axis(fig[2,1])
ax4 = Axis(fig[2,2])

case = collect(keys(out))[1]
caseout = out[case]
successful = caseout["error"] .< 0.01
h = hist!(ax1, caseout["error"], bins=0:0.01:1.5, label=case, color=(colorlist[6],1))
hist!(ax1, caseout["error"][successful], bins=h.bins, label=case, color=(colorlist[3],1))

h = hist!(ax2, caseout["treesize"], bins=0:500:13000, color=(colorlist[6],1))
hist!(ax2, caseout["treesize"][successful], bins=h.bins, color=(colorlist[3],1))

h = hist!(ax3, caseout["time"], bins=0:20:660, color=(colorlist[6],1))
hist!(ax3, caseout["time"][successful], bins=h.bins, color=(colorlist[3],1))

h = hist!(ax4, caseout["cost"], bins=0:40:1400, color=(colorlist[6],1))
hist!(ax4, caseout["cost"][successful], bins=h.bins, color=(colorlist[3],1))

ax1.xlabel = L"Final $d$"
ax1.ylabel = ax2.ylabel = ax3.ylabel = ax4.ylabel = "# Runs"
ax2.xlabel = "# Nodes"
ax3.xlabel = "Execution Time [s]"
ax4.xlabel = "Maneuver Length [s]"

save("output/bonsaipd_benchmark.pdf", fig, pt_per_unit=1/2.5)
