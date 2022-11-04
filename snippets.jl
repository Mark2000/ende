
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
