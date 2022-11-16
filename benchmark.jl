include("ende.jl")
include("endeplots.jl")

using FileIO
using Polynomials

N = 100

cases = Dict{String,Dict{Symbol,Any}}()
# cases["basicrrt"] = Dict(:pbend => 0., :biasfactor => Inf)
# cases["biasrrt"] = Dict(:pbend => 0.)
# cases["bonsairrt"] = Dict(:n => 20000)
# cases["thickbonsairrt"] = Dict(:globalbranchimprovement => false, :n => 15000)
cases["bonsaipdrrt"] = Dict(:samplecontrol => ()->[rand()*0.005; rand()*0.03], :ϵ => 0.01)
#:controlfnfactory => planarpdcontrolfactory,
function mccase()
    case = Dict{Symbol,Any}()
    # α = rand()*100
    # case[:α] = α
    # case[:distfn] = (x1, x2)->dist(x1, x2, α)
    # tmax = rand()*50
    # case[:tmax] = tmax
    # case[:sampletime] = ()->rand()*tmax
    # case[:pbias] = rand()
    # case[:kbias] = 25
    # case[:biasfactor] = sample([rand()*100, Inf], Weights([3/4,1/4]))

    # case[:pbend] = rand()*0.1
    # case[:kbend] = 50
    # case[:bendfactor] = sample([rand()*100, Inf], Weights([3/4,1/4]))
    case[:tmod] = sample([rand()*0.5, nothing], Weights([2/3,1/3]))
    case[:kmod] = sample([rand()*0.5, nothing], Weights([1/3,2/3]))
    # case[:globalbranchimprovement] = sample([true, false], Weights([2/3,1/3]))
    return case
end

Nmc = 200
for i in 1:Nmc
    cases[string(i)] = mccase()
end

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
    m_max = 0.66 # A*m^2
    J = [41.87 0.0 0.0; 0.0 41.87 0.0; 0.0 0.0 6.67] * 1e-3 # kg*m^2

    ωmax = 5*π/180 # rad/s

    # Alg. Params
    n = 10000
    ϵ = 0.0
    controlfnfactory = planarpdcontrolfactory
    samplestate() = ComponentArray(q = toscalarlast(randquat()), ω = randinunit(3) * ωmax)
    samplecontrol() = randinunit(2)
    sampletime() = rand()*25
    distfn(x1, x2) = dist(x1, x2, 30)

    pbias = 0.66
    kbias = 25
    biasfactor = 20

    pbend = 0.03
    kbend = 50
    bendfactor = 20
    tmod = 0.1
    kmod = nothing
    globalbranchimprovement = true

    for (sym, val) in values
        if isa(val,Function)
            @eval $(sym)(x...) = val(x...)
        else
            @eval $sym = $val
        end
    end

    Threads.@threads for i in 1:N
        # prob = AttitudeProblem(Bmag, m_max, J, ωmax, randquat(), Keepout[])
        prob = AttitudeProblem(Bmag.*[0.,0.,1.], m_max, J, ωmax, Quaternion(0.,0.,0.,1.), [Keepout([0.,0.,-1.],[-1.,0.,0.],50*π/180),Keepout(normalize([1.,1.,0.]),[0.,-1.,0.],20*π/180)])
        t0 = time()
        path, roadmap, states, controls = bonsairrt(prob, controlfnfactory, samplecontrol, samplestate, sampletime, distfn, n, ϵ, pbias, kbias, biasfactor, pbend, kbend, bendfactor; tmod=tmod, kmod=kmod, globalbranchimprovement=globalbranchimprovement)
        telapsed = time()-t0

        for (output, expr) in outputs
            if output == "time"
                out[case][output][i] = telapsed # TODO TODO TODO FIX THIS
            else
                out[case][output][i] = expr(path, roadmap, states, controls, prob, distfn)
            end
        end
    end
end
# FileIO.save("mcsweep_branchmod.jld2","run",[out,cases])

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


# Case Compare
fig = Figure(resolution=(800,800).*2, fontsize=24)
# fig = Figure(resolution=(1280,1280), fontsize=12)
ax1 = Axis(fig[1,1])
ax2 = Axis(fig[1,2])
ax3 = Axis(fig[2,1])
ax4 = Axis(fig[2,2])

exclude = ["biasrrt","basicrrt", "thickbonsairrt"]
colors = ax1.palette.patchcolor[]
for (i, (case, caseout)) in enumerate(out)
    if case ∉ exclude
        h = hist!(ax1, caseout["error"],bins=(0:0.01:1.5), normalization=:probability, label=case, color=colors[i])
        hist!(ax2, caseout["treesize"],bins=(13500:125:20000), normalization=:probability, color=colors[i])
        hist!(ax3, caseout["time"],bins=(0:15:500), normalization=:probability, color=colors[i])
        hist!(ax4, caseout["cost"],bins=(0:30:700), normalization=:probability, color=colors[i])
    end
end

ax1.xlabel = "Final Error"
ax1.ylabel = ax2.ylabel = ax3.ylabel = "% Runs"
ax2.xlabel = "Tree Nodes"
ax3.xlabel = "Execution Time [s]"
ax4.xlabel = "Solution Cost [s]"

fig[1,3] = Legend(fig, ax1)
save("output/comp_bonsaipd.png", fig)
