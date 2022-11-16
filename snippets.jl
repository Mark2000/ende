
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


    # qdot_0 = (1/2* Ξ(q_lin) * ω_0)[1:3]
    # ωdot_0 = - prob.Jinv * ω_0 × (prob.J * ω_0)
    #
    # J1 = prob.J[1,1]
    # J2 = prob.J[2,2]
    # J3 = prob.J[3,3]
    # dwdw = [0 (J2-J3)/J1*ω_0[3] (J2-J3)/J1*ω_0[2];
    #         (J3-J1)/J2*ω_0[3] 0 (J3-J1)/J2*ω_0[1]
    #         (J1-J2)/J3*ω_0[2] (J1-J2)/J3*ω_0[1] 0];
    #
    # show(stdout, "text/plain", dwdw)
    #
    # dqdy = [0 ω_0[3] -ω_0[2] q_lin.s -q_lin.v3 q_lin.v2;
    #         -ω_0[3] 0 ω_0[1] q_lin.v3 q_lin.s -q_lin.v1;
    #         ω_0[2] -ω_0[1] 0 -q_lin.v2 q_lin.v1 q_lin.s]/2
    #
    # show(stdout, "text/plain", dqdy)
    #
    # A = [dqdy qdot_0;
    #      zeros(3,3) dwdw ωdot_0;
    #      zeros(1,7)]
    #
    # rotm = QuatRotation(q_BI0)
    # B_B = rotm * prob.B(tinit)
    # xI_B = rotm * [1.0, 0.0, 0.0]
    # Cx_B = normalize((B_B × xI_B) × B_B)
    # Cy_B = normalize(B_B × Cx_B)

    # B = [zeros(3,2);
    #      B_B×Cx_B./diag(J) B_B×Cy_B./diag(J);
    #      zeros(1,2)]

    # show(stdout, "text/plain", A)
    # show(stdout, "text/plain", B)
    #


    function lqrcontrolfactory(prob::AttitudeProblem, tinit::Real, xinit::State, xtarget::State)
        # state vector y = [qx, qy, qz, wx, wy, wz, 1]
        q_FI = normalize(fromscalarlast(xtarget.q))
        q_BI0 = normalize(fromscalarlast(xinit.q))
        q_BF0 = q_BI0*inv(q_FI) # initial q error
        q_lin = q_BF0

        ω_0 = xinit.ω


        A = zeros(6,6)
        A[1:3,4:6] = 0.5I(3)

        B = [I(3); prob.Jinv]

        Q = 1.0I(6)
        Q[4,4] = Q[5,5] = Q[6,6] = 1
        R = 0.001I(3)

        L = lqr(Continuous,A,B,Q,R)#*1e-6

        # @info "LQR out"
        @show L

        function lqrcontrol(t::Real, k::Vector, x::State)::Vector
            q_BI = normalize(fromscalarlast(x.q))
            rotm = QuatRotation(q_BI)
            B_B = rotm * B_I

            q_BF = q_BI*inv(q_FI)
            y = [q_BF.v1; q_BF.v2; q_BF.v3; x.ω] # TODO omega wrong

            T = -L*y
            Tonplane = T - T⋅B_B./norm(B_B)*B_B./norm(B_B)

            m_mag = min(norm(T)/norm(B_B),m_max)
            # @show m_mag m_max

            m_B = m_mag * (B_B./norm(B_B)) × (Tonplane./norm(Tonplane))
            return m_B × B_B
        end
        return lqrcontrol
    end
