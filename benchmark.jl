include("ende.jl")
include("endeplots.jl")

using FileIO
using Dates
# using Polynomials

N = 8

cases = Dict{String,Dict{Symbol,Any}}()
timeout = 200

cases["rrt-const"] = Dict(:controlfnfact => planarcontrolfactory, :samplecontrol => ()->randinunit(2), :timeout => timeout, :pbend => 0, :biasfactor => Inf)
cases["rrt-pd"] = Dict(:controlfnfact => planarpdcontrolfactory, :samplecontrol => ()->[rand()*0.005; rand()*0.03], :timeout => timeout, :pbend => 0, :biasfactor => Inf)
#
# cases["softmax-const"] = Dict(:controlfnfact => planarcontrolfactory, :samplecontrol => ()->randinunit(2), :timeout => timeout, :pbend => 0)
# cases["softmax-pd"] = Dict(:controlfnfact => planarpdcontrolfactory, :samplecontrol => ()->[rand()*0.005; rand()*0.03], :timeout => timeout, :pbend => 0)

cases["bonsai-const"] = Dict(:controlfnfact => planarcontrolfactory, :samplecontrol => ()->randinunit(2), :timeout => timeout)
cases["bonsai-pd"] = Dict(:controlfnfact => planarpdcontrolfactory, :samplecontrol => ()->[rand()*0.005; rand()*0.03], :timeout => timeout)

# cases["bonsaibiasmax-pd"] = Dict(:controlfnfact => planarpdcontrolfactory, :samplecontrol => ()->[rand()*0.005; rand()*0.03], :timeout => timeout, :biasfactor => Inf)
# cases["bonsaibendmax-pd"] = Dict(:controlfnfact => planarpdcontrolfactory, :samplecontrol => ()->[rand()*0.005; rand()*0.03], :timeout => timeout, :bendfactor => Inf)
#
# cases["bonsaibiasbendmax-pd"] = Dict(:controlfnfact => planarpdcontrolfactory, :samplecontrol => ()->[rand()*0.005; rand()*0.03], :timeout => timeout, :biasfactor => Inf, :bendfactor => Inf)


outputs = Dict{String,Union{Symbol,Function}}()
outputs["time"] = ()->nothing
outputs["cost"] = (path, roadmap, states, controls, prob, distfn)->timecost(controls, path)
outputs["qerror"] = (path, roadmap, states, controls, prob, distfn)->dist(states[path[end]], stategoal(prob), 0)
outputs["werror"] = (path, roadmap, states, controls, prob, distfn)->dist(states[path[end]], stategoal(prob), 1e10)/1e10
outputs["error"] = (path, roadmap, states, controls, prob, distfn)->dist(states[path[end]], stategoal(prob), 30)
outputs["treesize"] = (path, roadmap, states, controls, prob, distfn)->length(states)

out = Dict{String,Dict{String,Vector}}(case => Dict(output => zeros(N) for output in keys(outputs)) for case in keys(cases))

for case in collect(keys(cases))
    values = cases[case]
    @info "Benchmarking" case values
    # Sys. Params
    Bmag = 5e-5 # Teslas (kg/As^2)
    B_hat = [0.; 0.; 1.]
    B_I = B_hat * Bmag

    m_max = 0.66 # A*m^2
    J = [41.87 0.0 0.0; 0.0 41.87 0.0; 0.0 0.0 6.67] * 1e-3 # kg*m^2

    ωmax = 5*π/180 # rad/s

    # Alg. Params
    n = 1000000 # using timeout
    ϵ = 0.0
    samplestate() = ComponentArray(q = toscalarlast(randquat()), ω = randinunit(3) * ωmax)
    samplecontrol() = randinunit(2)
    sampletime() = rand()*25
    distfn(x1, x2) = dist(x1, x2, 30)
    controlfnfact(x...) = planarpdcontrolfactory(x...)

    pbias = 0.66
    kbias = 25
    biasfactor = 20

    pbend = 0.03
    kbend = 50
    bendfactor = 20
    tmod = 0.1
    kmod = nothing
    globalbranchimprovement = true

    timeout = Inf

    for (sym, val) in values
        if isa(val,Function)
            @eval $(sym)(x...) = $(val)(x...)
        else
            @eval $sym = $val
        end
    end

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
# groupings = [["rrt-const","rrt-pd"],["softmax-const","softmax-pd"],["bonsai-const","bonsai-pd"]]
groupings = [["bonsai-pd"],["bonsaibiasmax-pd"],["bonsaibendmax-pd"],["bonsaibiasbendmax-pd"]]
results = ["error","cost","treesize"]
fig = Figure()
for (axi, result) in enumerate(results)
    ax = Axis(fig[1,axi])
    xs = vcat([repeat([g],length(out[trial][result])) for (g,group) in enumerate(groupings) for (t,trial) in enumerate(group) ]...)
    ys = vcat([out[trial][result] for (g,group) in enumerate(groupings) for (t,trial) in enumerate(group) ]...)
    dodge = vcat([repeat([t],length(out[trial][result])) for (g,group) in enumerate(groupings) for (t,trial) in enumerate(group) ]...)
    boxplot!(ax, xs, ys, dodge = dodge, color = dodge) # show_notch = true
    ax.xlabel = "Trial"
    ax.ylabel = result
end



# Case Compare
fig = Figure(resolution=(800,800).*2, fontsize=24)
# fig = Figure(resolution=(1280,1280), fontsize=12)
ax1 = Axis(fig[1,1])
ax2 = Axis(fig[1,2])
ax3 = Axis(fig[2,1])
ax4 = Axis(fig[2,2])

exclude = []
colors = ax1.palette.patchcolor[]
for (i, (case, caseout)) in enumerate(out)
    if case ∉ exclude
        h = hist!(ax1, caseout["error"], normalization=:probability, label=case, color=colors[i])
        hist!(ax2, caseout["treesize"], normalization=:probability, color=colors[i])
        hist!(ax3, caseout["time"], normalization=:probability, color=colors[i])
        hist!(ax4, caseout["cost"], normalization=:probability, color=colors[i])
    end
end

ax1.xlabel = "Final Error"
ax1.ylabel = ax2.ylabel = ax3.ylabel = "% Runs"
ax2.xlabel = "Tree Nodes"
ax3.xlabel = "Execution Time [s]"
ax4.xlabel = "Solution Cost [s]"

fig[1,3] = Legend(fig, ax1)
save("output/comp_bonsaipd.png", fig)


#
# function mccase()
#     case = Dict{Symbol,Any}()
#     # α = rand()*100
#     # case[:α] = α
#     # case[:distfn] = (x1, x2)->dist(x1, x2, α)
#     # tmax = rand()*50
#     # case[:tmax] = tmax
#     # case[:sampletime] = ()->rand()*tmax
#     # case[:pbias] = rand()
#     # case[:kbias] = 25
#     # case[:biasfactor] = sample([rand()*100, Inf], Weights([3/4,1/4]))
#
#     # case[:pbend] = rand()*0.1
#     # case[:kbend] = 50
#     # case[:bendfactor] = sample([rand()*100, Inf], Weights([3/4,1/4]))
#     case[:tmod] = sample([rand()*0.5, nothing], Weights([2/3,1/3]))
#     case[:kmod] = sample([rand()*0.5, nothing], Weights([1/3,2/3]))
#     # case[:globalbranchimprovement] = sample([true, false], Weights([2/3,1/3]))
#     return case
# end
#
# Nmc = 200
# for i in 1:Nmc
#     cases[string(i)] = mccase()
# end


# MC Compare: Lumped
fig = Figure()
ax1 = Axis(fig[1,1])
ax2 = Axis(fig[1,2])
ax3 = Axis(fig[2,1])
ax4 = Axis(fig[2,2])
ax5 = Axis(fig[3,1])
ax6 = Axis(fig[3,2])

hist!(ax1, vcat([case["error"] for case in values(out)]...), normalization=:probability)
ax1.xlabel = "Final Error"

hist!(ax2, vcat([case["treesize"] for case in values(out)]...), normalization=:probability)
ax2.xlabel = "Tree Nodes"

hist!(ax3, vcat([case["qerror"]*180/pi for case in values(out)]...), normalization=:probability)
ax3.xlabel = "Angle Error [deg]"

hist!(ax4, vcat([case["werror"]*180/pi for case in values(out)]...), normalization=:probability)
ax4.xlabel = "Angle Rate Error [deg]"

hist!(ax5, vcat([case["time"] for case in values(out)]...), normalization=:probability)
ax5.xlabel = "Execution Time [s]"

hist!(ax6, vcat([case["cost"] for case in values(out)]...), normalization=:probability)
ax6.xlabel = "Solution Cost [s]"

ax1.ylabel = ax2.ylabel = ax3.ylabel = ax4.ylabel = ax5.ylabel = ax6.ylabel = "% Runs"

# MC Compare: Correlation
fig = Figure()
cols = 4
errorout = vcat([mean(case["error"]) for case in values(out)]...)
derrorout = vcat([std(case["error"]) for case in values(out)]...)*1.0
for (i, sym) in enumerate([sym for (sym,val) in cases["1"] if isa(val,Union{Real,Bool,Nothing})])
    ax = Axis(fig[Int(ceil(i/cols)), mod1(i,cols)])
    ax.xlabel = string(sym)

    vals = Any[params[sym] for (case, params) in cases]
    numeric = (.!isnothing.(vals)) .& .!(vals .== Inf)
    # vals[isnothing.(vals)] .= NaN
    # vals[isinf.(vals)] .= NaN
    valnum = convert(Array{Float64,1}, vals[numeric])
    scatter!(ax, valnum, errorout[numeric])
    errorbars!(ax, valnum, errorout[numeric], derrorout[numeric], color=:grey)

    order = 2
    p = Polynomials.fit(valnum, errorout[numeric], order)
    lines!(ax, sort(valnum), p.(sort(valnum)))

    p = Polynomials.fit(valnum, errorout[numeric]+derrorout[numeric], order)
    lines!(ax, sort(valnum), p.(sort(valnum)))

    p = Polynomials.fit(valnum, errorout[numeric]-derrorout[numeric], order)
    lines!(ax, sort(valnum), p.(sort(valnum)))

    hlines!(ax, [mean(errorout[.!numeric])], color = :red)
end
