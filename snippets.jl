
#
# n = 5000
# ϵ = 0.10
# # pbiases = [0.05, 0.1, 0.25]
# pbias = 0.10
# dt = 10
# # ks = [5, 25, 100]
# k = 25
# α = 5
# f1() = randinunit(2)
# f2() = randunit(2)
# samplecontrols = [f1, f2]
#
# ntrials = 50
# lengths = Dict()
# errors = Dict()
# times = Dict()
# for i in 1:2
#     lengths[i] = zeros(ntrials)
#     errors[i] = zeros(ntrials)
#     times[i] = zeros(ntrials)
# end
#
# Threads.@threads for i in 1:ntrials
#     prob = AttitudeProblem(B_I, m_max, J)
#     for (f, samplecontrol) in enumerate(samplecontrols)
#         @info "Testing" samplecontrol i
#         t0 = time()
#         path, roadmap, states, controls = dynamicrrt(prob, samplecontrol, DIST_TODO samplestate, n, ϵ, pbias, dt, k)
#         telapsed = time()-t0
#         lengths[f][i] = length(path)
#         errors[f][i] = dist(states[path[end]], stategoal(prob))
#         times[f][i] = telapsed
#     end
# end
#
# fig = Figure()
# ax1 = Axis(fig[1,1])
# xs = collect(Iterators.flatten([repeat([k*1.0], length(lengths[k])) for k in 1:2]))
# ys = collect(Iterators.flatten([lengths[k]*1.0 for k in 1:2]))
# boxplot!(ax1, xs, ys, width=0.5)
# ax1.xlabel = "f"
# ax1.ylabel = "length"
# # ax1.xticks = pbiases
#
# ax2 = Axis(fig[2,1])
# xs = collect(Iterators.flatten([repeat([k*1.0], length(errors[k])) for k in 1:2]))
# ys = collect(Iterators.flatten([errors[k]*1.0 for k in 1:2]))
# boxplot!(ax2, xs, ys, width=0.5)
# ax2.xlabel = "f"
# ax2.ylabel = "error"
# # ax2.xticks = pbiases
#
# ax3 = Axis(fig[3,1])
# xs = collect(Iterators.flatten([repeat([k*1.0], length(times[k])) for k in 1:2]))
# ys = collect(Iterators.flatten([times[k]*1.0 for k in 1:2]))
# boxplot!(ax3, xs, ys, width=0.5)
# ax3.xlabel = "f"
# ax3.ylabel = "time"
# # ax3.xticks = pbiases


function sst(
    prob,
    samplecontrol::Function,
    samplestate::Function,
    sampletime::Function,
    dist::Function,
    n,
    ϵ,
    δbn,
    δₛ,
    pbias,
)
    # Currently assumes time is cost
    xstart = statestart(prob)
    xgoal = stategoal(prob)

    states = [xstart]
    controls = [ControlCommand(0.0, [0.0; 0.0])]
    path = []
    roadmap = Dict()

    Vactive = [1]
    Vinactive = Integer[]
    S = Dict{Int,Union{Int,Nothing}}(1 => 1)

    function cost(istate)
        path = pathto(roadmap, istate)
        sum([control.t for control in controls[path]])
    end

    imin = 1
    for iter = 1:n
        if mod(iter, 1000) == 0
            @info "on step" iter
        end

        # Best First SST
        if rand() < pbias
            xrand = xgoal
        else
            xrand = samplestate()
        end
        Inear = [i for i in Vactive if dist(xrand, states[i]) < δbn]
        isel =
            isempty(Inear) ?
            Vactive[argmin([dist(x, xrand) for x in states[Vactive]])] :
            Inear[argmin([cost(i) for i in Inear])]
        xsel = states[isel]

        u = samplecontrol()
        dt = sampletime()
        xnew, valid = stepstate(xsel, u, dt, prob)

        if valid
            push!(states, xnew)
            push!(controls, ControlCommand(dt, u))
            ixnew = length(states)

            # Is node locally the best
            isnodelocalbest = false
            isnew = collect(keys(S))[argmin([
                dist(x, xnew) for x in states[collect(keys(S))]
            ])]
            if dist(states[isnew], xnew) > δₛ
                isnew = ixnew
                S[isnew] = nothing
            end
            ixpeer = S[isnew]
            if isnothing(ixpeer) || cost(isel) + dt < cost(ixpeer)
                isnodelocalbest = true
            end

            if isnodelocalbest
                push!(Vactive, ixnew)
                roadmap[ixnew] = isel

                if dist(xnew, xgoal) < dist(states[imin], xgoal) ||
                   ((dist(xnew, xgoal) < ϵ) && (cost(ixnew) < cost(imin)))
                    imin = length(states)
                    if xgoal == xrand
                        @info "Biased step"
                    end
                    @info xnew imin dist(xnew, xgoal) cost(imin)
                end

                # Prune dominated nodes
                if !isnothing(ixpeer)
                    filter!(x -> x ≠ ixpeer, Vactive)
                    push!(Vinactive, ixpeer)
                end
                S[isnew] = ixnew
                while !isnothing(ixpeer) &&
                          ixpeer ∉ values(roadmap) &&
                          ixpeer ∈ Vinactive
                    ixparent = roadmap[ixpeer]
                    delete!(roadmap, ixpeer)
                    filter!(x -> x ≠ ixpeer, Vinactive)
                    ixpeer = ixparent
                end
            end
        end
    end

    # if dist(xnew, xgoal) < ϵ
    # @warn "Did not converge within ϵ"
    path = pathto(roadmap, imin)

    @info "Path length:" sum([control.t for control in controls[path]])

    return path, roadmap, states, controls
end
