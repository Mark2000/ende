using ComponentArrays
using DifferentialEquations
using LinearAlgebra
using Optim
using Quaternions
using Random
using Rotations
using StatsBase

Quaternion = Quaternions.Quaternion

# Random Element Generation Functions
randquat() = normalize(Quaternion(randn(4)...)) # random unit quaternion

randunit(n) = normalize([randn(n)...]) # random point on n-sphere

function randinunit(n) # random point in n-sphere
    while true
        v = 2 * rand(n) .- 1
        if norm(v) ≤ 1
            return v
        end
    end
end

# Utilities
fromscalarlast(q) = Quaternion(q[4], q[1], q[2], q[3])
toscalarlast(q) = [q.v1; q.v2; q.v3; q.s]

Ξ(q) = [q.s -q.v3 q.v2; q.v3 q.s -q.v1; -q.v2 q.v1 q.s; -q.v1 -q.v2 -q.v3]

function angle(a, b)
    return acos(clamp(a ⋅ b / (norm(a) * norm(b)), -1, 1))
end

softmax(x) = exp.(x) ./ sum(exp.(x))

# Distance Metrics
function dist(x1, x2, α = 5)
    """
    Distance metric between (q, ω) states
        eigangle(q₁,q₂) + α * norm(ω₁-ω₂)
    """

    qw = max(min(x1.q ⋅ x2.q, 1.0), -1.0)
    dangle = abs(2 * acos(qw))

    dω = norm(x1.ω - x2.ω)

    dangle + α * dω
end

function distB(x1, x2, B, α = 5, ϵ = 1, β = 1)
    """
    Distance metric between (q, ω) states that penalizes out of B plane error
        eigangle(q₁,q₂) * (1 + ϵ * eigangleₒₒₚ(q₁,q₂))
        + α * norm(ω₁-ω₂) * (1 + β * normₒₒₚ(ω₁-ω₂))
    """

    if ϵ == β == 0
        return dist(x1, x2, α)
    end

    B_hat = normalize(B)
    qw = max(min(x1.q ⋅ x2.q, 1.0), -1.0)
    dangle = abs(2 * acos(qw)) # in range [0,π]
    qv = -x1.q[1:3] .* x2.q[4] + x1.q[4] .* x2.q[1:3] - x1.q[1:3] × x2.q[1:3]
    eigax = qv / sqrt(1 - qw^2)
    oopangle = abs(B_hat ⋅ eigax) # in range [0,1]

    dω = norm(x1.ω - x2.ω) # in range [0,ωmax]
    oopω = abs(B_hat ⋅ normalize(x1.ω - x2.ω))  # in range [0,1]
    if isnan(oopω)
        oopω = 0
    end

    cost = dangle * (1 + ϵ * oopangle) + α * dω * (1 + β * oopω)

    # @info cost dangle oopangle dω oopω
    return cost
end

# Problem Definition Elements
struct Keepout
    sensor::Vector
    obstacle::Vector
    halfangle::Real
end

struct AttitudeProblem
    B::Vector
    m_max::Real
    J::Matrix
    Jinv::Matrix
    qstart::Quaternion
    ωstart::Vector
    qgoal::Quaternion
    ωgoal::Vector
    ωmax::Real
    keepouts::Vector{Keepout}
end

AttitudeProblem(
    B::Vector,
    m_max::Real,
    J::Matrix,
    qstart::Quaternion,
    ωstart::Vector,
    qgoal::Quaternion,
    ωgoal::Vector,
    ωmax::Real,
    keepouts::Vector{Keepout},
) = AttitudeProblem(
    B,
    m_max,
    J,
    inv(J),
    qstart,
    ωstart,
    qgoal,
    ωgoal,
    ωmax,
    keepouts,
)

AttitudeProblem(
    B::Real,
    m_max::Real,
    J::Matrix,
    qstart::Quaternion,
    ωstart::Vector,
    qgoal::Quaternion,
    ωgoal::Vector,
    ωmax::Real,
    keepouts::Vector{Keepout},
) = AttitudeProblem(
    B * randunit(3),
    m_max,
    J,
    qstart,
    ωstart,
    qgoal,
    ωgoal,
    ωmax,
    keepouts,
)

AttitudeProblem(B, m_max::Real, J::Matrix, ωmax::Real) = AttitudeProblem(
    B,
    m_max,
    J,
    randquat(),
    [0.0; 0.0; 0.0],
    Quaternion(1.0),
    [0.0; 0.0; 0.0],
    ωmax,
    Keepout[],
)
AttitudeProblem(B, m_max::Real, J::Matrix, ωmax::Real, qstart::Quaternion) =
    AttitudeProblem(
        B,
        m_max,
        J,
        qstart,
        [0.0; 0.0; 0.0],
        Quaternion(1.0),
        [0.0; 0.0; 0.0],
        ωmax,
        Keepout[],
    )
AttitudeProblem(
    B,
    m_max::Real,
    J::Matrix,
    ωmax::Real,
    qstart::Quaternion,
    keepouts::Vector{Keepout},
) = AttitudeProblem(
    B,
    m_max,
    J,
    qstart,
    [0.0; 0.0; 0.0],
    Quaternion(1.0),
    [0.0; 0.0; 0.0],
    ωmax,
    keepouts,
)


statestart(prob::AttitudeProblem) =
    ComponentArray(q = toscalarlast(prob.qstart), ω = prob.ωstart)
stategoal(prob::AttitudeProblem) =
    ComponentArray(q = toscalarlast(prob.qgoal), ω = prob.ωgoal)

struct ControlledAttitudeProblem
    prob::AttitudeProblem
    u::Function
end

struct ControlCommand
    t::Real
    u::Vector
end

# System Dynamics
function eulereqns(dx, x, params, t)
    """
    Derivative of Euler equations with intertially fixed dipole in constant B field
    """
    J = params.prob.J
    Jinv = params.prob.Jinv
    B_I = params.prob.B
    q = normalize(fromscalarlast(x.q))
    rotm = QuatRotation(q)
    B_B = rotm * B_I

    xI_B = rotm * [1.0, 0.0, 0.0]
    Bx_B = normalize((B_B × xI_B) × B_B)
    By_B = normalize(B_B × Bx_B)
    m_B = params.prob.m_max * (params.u(t)[1] * Bx_B + params.u(t)[2] * By_B)

    L_B = m_B × B_B

    α = Jinv * (L_B - x.ω × (J * x.ω))
    qdot = 1 / 2 * Ξ(q) * x.ω

    dx.q = qdot
    dx.ω = α
    return nothing
end

function stepstate(state, u, t, prob::AttitudeProblem)
    ufun(t) = u
    params = ControlledAttitudeProblem(prob, ufun)
    ode = ODEProblem(eulereqns, state, (0, t), params)
    sol = solve(ode, reltol = 1e-7, abstol = 1e-9) # save_everystep = false
    for x in sol.u
        if norm(x.ω) > prob.ωmax
            return sol.u[end], false
        end
        if length(prob.keepouts) > 0
            rotm = QuatRotation(normalize(fromscalarlast(x.q)))
            for keepout in prob.keepouts
                obs_B = rotm * keepout.obstacle
                if angle(keepout.sensor, obs_B) < keepout.halfangle
                    return sol.u[end], false
                end
            end
        end
    end
    sol.u[end], true
end

# Roadmap Traversal
function pathto(roadmap, istate::Integer)
    path = [istate]
    while path[1] != 1
        pushfirst!(path, roadmap[path[1]])
    end
    return path
end

timecost(roadmap, controls, istate::Integer) =
    sum([control.t for control in controls[pathto(roadmap, istate)]])

timecost(controls, path) = sum([control.t for control in controls[path]])

# Search Functions
function dynamicrrt(
    prob,
    samplecontrol::Function,
    samplestate::Function,
    sampletime::Function,
    dist::Function,
    n,
    ϵ,
    pbias,
    k,
)
    xstart = statestart(prob)
    xgoal = stategoal(prob)

    states = [xstart]
    controls = [ControlCommand(0.0, [0.0; 0.0])]
    path = []
    roadmap = Dict()
    imin = 1
    for iter = 1:n
        if mod(iter, 1000) == 0
            @info "on step" iter
        end

        unew = [0.0, 0.0]
        xnew = ComponentArray()
        valid = false
        inear = imin
        dt = 0.0

        if rand() < pbias
            xrand = xgoal
            xnear = states[imin]
            dcurrent = Inf
            for i = 1:k
                utest = samplecontrol()
                dtcurrent = sampletime()
                xtest, valid = stepstate(xnear, utest, dtcurrent, prob)
                dtest = dist(xtest, xgoal)
                if valid & (dist(xtest, xgoal) < dcurrent)
                    xnew = xtest
                    unew = utest
                    dcurrent = dtest
                    dt = dtcurrent
                end
            end
        else
            xrand = samplestate()
            inear = argmin([dist(x, xrand) for x in states])
            xnear = states[inear]

            unew = samplecontrol()
            dt = sampletime()
            xnew, valid = stepstate(xnear, unew, dt, prob)
        end

        if valid
            push!(states, xnew)
            push!(controls, ControlCommand(dt, unew))
            roadmap[length(states)] = inear
            if dist(xnew, xgoal) < dist(states[imin], xgoal)
                imin = length(states)
                if xrand == xgoal
                    @info "Biased step found new minimum:"
                end
                @info xnew imin dist(xnew, xgoal)
            end
            if dist(xnew, xgoal) < ϵ
                @info "Converged within ϵ"
                path = pathto(roadmap, length(states))
                break
            end
        end
    end

    if path == []
        @warn "Did not converge within ϵ"
        path = pathto(roadmap, imin)
    end

    @info "Path time:" sum([control.t for control in controls[path]])

    return path, roadmap, states, controls
end

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

function bonsairrt(
    prob,
    samplecontrol::Function,
    samplestate::Function,
    sampletime::Function,
    dist::Function,
    n,
    ϵ,
    pbias,
    kbias,
    biasfactor,
    pbend,
    kbend,
    bendfactor;
    tmod = nothing,
    umod = nothing,
    globalbranchimprovement = false,
)
    xstart = statestart(prob)
    xgoal = stategoal(prob)

    states = [xstart]
    controls = [ControlCommand(0.0, [0.0; 0.0])]
    path = []
    roadmap = Dict()
    imin = 1 # Index of best state

    for iter = 1:n
        if mod(iter, 1000) == 0
            @info "on step" iter
        end

        if imin > 1 && rand() < pbend
            # Choose to bend the best branch if bendfactor = Inf, otherwise pick a random branch
            inear = imin
            if bendfactor != Inf
                leaves = setdiff(Set(keys(roadmap)), Set(values(roadmap)))
                while true # Pick a random leave based on distribution
                    weight =
                        softmax(-[dist(x, xgoal) for x in states] * bendfactor)
                    inear = sample(1:length(states), Weights(weight))
                    if inear in leaves # must be a leaf so we test a whole branch
                        break
                    end
                end
            end
            branch = pathto(roadmap, inear)

            # Set the threshold for saving a new branch
            dcurrent = nothing
            if globalbranchimprovement
                dcurrent = dist(states[imin], xgoal)
            else
                dcurrent = min([dist(x, xgoal) for x in states[branch]]...)
            end

            # Test at most kbranch branches
            dminbranch = nothing
            beststates = nothing
            bestcontrols = nothing
            bestelbow = nothing
            kmin = nothing
            jelbow = nothing
            for iter = 1:kbend
                jelbow = rand(2:length(branch))
                branchcontrols = controls[branch[jelbow+1:end]]

                tnew = sampletime()
                unew = samplecontrol()
                if !isnothing(umod)
                    unew = controls[branch[jelbow]].u .+ umod .* randinunit(2)
                    unew = norm(unew) > 1 ? unew / norm(unew) : unew
                end
                if !isnothing(tmod)
                    tnew =
                        controls[branch[jelbow]].t *
                        (1 - tmod + 2 * rand() * tmod)
                end
                controltest = ControlCommand(tnew, unew)

                xfirst, valid = stepstate(
                    states[branch[jelbow-1]],
                    controltest.u,
                    controltest.t,
                    prob,
                )
                if valid
                    branchstates = [xfirst]
                    for control in branchcontrols
                        xnext, valid = stepstate(
                            branchstates[end],
                            control.u,
                            control.t,
                            prob,
                        )
                        if valid
                            push!(branchstates, xnext)
                        else
                            break
                        end
                    end
                    pushfirst!(branchcontrols, controltest)
                    #TODO delete old branch?
                    branchdist = [dist(x, xgoal) for x in branchstates]
                    dminbranch = min(branchdist...)
                    if dminbranch < dcurrent
                        @info "control update" controls[branch[jelbow]] controltest
                        beststates = branchstates
                        bestcontrols = branchcontrols
                        bestelbow = jelbow
                        dcurrent = dminbranch
                        kmin = argmin(branchdist)
                        break
                    end
                end
            end
            if !isnothing(beststates)
                if dminbranch < dist(states[imin], xgoal)
                    @info "Bonsai reduced error" (jelbow, length(branch)) dist(
                        states[imin],
                        xgoal,
                    ) dcurrent beststates[kmin]
                    imin = length(states) + kmin
                end
                roadmap[length(states)+1] = branch[bestelbow-1]
                for i = length(states)+1:length(states)+length(beststates)-1
                    roadmap[i+1] = i
                end
                append!(states, beststates)
                append!(controls, bestcontrols)
            end

        else
            unew = [0.0, 0.0]
            xnew = ComponentArray()
            valid = false
            inear = imin
            dt = 0.0
            biased = false

            if rand() < pbias
                biased = true
                xrand = xgoal
                inear = imin
                if biasfactor != Inf
                    weight =
                        softmax(-[dist(x, xgoal) for x in states] * biasfactor)
                    inear = sample(1:length(states), Weights(weight))
                end
                xnear = states[inear]
                dcurrent = Inf
                for i = 1:kbias
                    utest = samplecontrol()
                    dtcurrent = sampletime()
                    xtest, valid = stepstate(xnear, utest, dtcurrent, prob)
                    dtest = dist(xtest, xgoal)
                    if valid & (dist(xtest, xgoal) < dcurrent) # TODO: how often to find new solution
                        xnew = xtest
                        unew = utest
                        dcurrent = dtest
                        dt = dtcurrent
                    end
                end
            else
                xrand = samplestate()
                inear = argmin([dist(x, xrand) for x in states])
                xnear = states[inear]

                unew = samplecontrol()
                dt = sampletime()
                xnew, valid = stepstate(xnear, unew, dt, prob)
            end

            if valid
                push!(states, xnew)
                push!(controls, ControlCommand(dt, unew))
                roadmap[length(states)] = inear
                if dist(xnew, xgoal) < dist(states[imin], xgoal) # todo move checks outside if
                    imin = length(states)
                    if biased
                        @info "Biased step found new minimum:"
                    else
                        @info "Random step found new minimum:"
                    end
                    @info xnew imin dist(xnew, xgoal)
                end
                if dist(xnew, xgoal) < ϵ
                    @info "Converged within ϵ"
                    path = pathto(roadmap, length(states))
                    break
                end
            end
        end
    end

    if path == []
        @warn "Did not converge within ϵ"
        path = pathto(roadmap, imin)
    end

    @info "Path time:" sum([control.t for control in controls[path]])

    return path, roadmap, states, controls
end

# Postprocessing Functions
function controlfnfactory(path, controls)
    tcontrols = cumsum([control.t for control in controls[path]])
    tmax = tcontrols[end]
    function ufun(t)
        if t == 0.0
            return controls[path[2]].u
        end
        icontrol = findfirst(t .< tcontrols)
        if isnothing(icontrol)
            return [0.0, 0.0]
        end
        return controls[path[icontrol]].u
    end
    return ufun, tmax, tcontrols
end

function resamplesolution(
    prob::AttitudeProblem,
    ufun::Function,
    tmax::Real;
    saveat::Real = 1,
    comparison = nothing,
    tstops = nothing,
)
    params = ControlledAttitudeProblem(prob, ufun)
    ode = ODEProblem(eulereqns, statestart(prob), (0, tmax), params)
    sol = solve(
        ode,
        reltol = 1e-7,
        abstol = 1e-9,
        saveat = saveat,
        tstops = tstops,
    )
    if !isnothing(comparison)
        @info "Integrated difference should be small" dist(
            sol.u[end],
            comparison,
        )
    end
    return sol
end
