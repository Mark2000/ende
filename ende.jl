using ComponentArrays: ComponentArray
using DifferentialEquations
using LinearAlgebra
using Optim
using Quaternions
using Random
using StatsBase
using ControlSystemsBase

Quaternion = Quaternions.Quaternion
function QuatRotation(q::Quaternion)
    qw = q.s
    qx = q.v1
    qy = q.v2
    qz = q.v3
    [1-2*qy^2-2*qz^2 2*qx*qy-2*qz*qw 2*qx*qz+2*qy*qw
     2*qx*qy+2*qz*qw 1-2*qx^2-2*qz^2 2*qy*qz-2*qx*qw
     2*qx*qz-2*qy*qw 2*qy*qz+2*qx*qw 1-2*qx^2-2*qy^2]
end

# Random Element Generation Functions
randquat() = Quaternion(normalize(randn(4))...) # random unit quaternion

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

State = typeof(ComponentArray(q = [0.,0.,0.,1.], ω = [0.,0.,0.]))

# Distance Metrics
function dist(x1::State, x2::State, α = 5)
    """
    Distance metric between (q, ω) states
        eigangle(q₁,q₂) + α * norm(ω₁-ω₂)
    """

    qw = max(min(x1.q ⋅ x2.q, 1.0), -1.0)
    dangle = abs(2 * acos(qw))

    dω = norm(x1.ω - x2.ω)

    dangle + α * dω
end

function distB(x1::State, x2::State, B, α = 5, ϵ = 1, β = 1)
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

# TODO: Base.@kwdef
struct AttitudeProblem
    B::Function
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
    B::Function,
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

function AttitudeProblem(
    B::Union{Real, Vector},
    m_max::Real,
    J::Matrix,
    qstart::Quaternion,
    ωstart::Vector,
    qgoal::Quaternion,
    ωgoal::Vector,
    ωmax::Real,
    keepouts::Vector{Keepout},
)
    if isa(B, Real)
        Bvec = B * randunit(3)
    else
        Bvec = B
    end
    AttitudeProblem(
        (t)->Bvec,
        m_max,
        J,
        qstart,
        ωstart,
        qgoal,
        ωgoal,
        ωmax,
        keepouts,
    )
end

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

statestart(prob::AttitudeProblem)::State =
    ComponentArray(q = toscalarlast(prob.qstart), ω = prob.ωstart)
stategoal(prob::AttitudeProblem)::State =
    ComponentArray(q = toscalarlast(prob.qgoal), ω = prob.ωgoal)

struct ControlCommand
    # Use controller u(t, k, x) for t seconds
    t::Real
    controlfn::Function
    k
end

ControlCommand() = ControlCommand(0.0, (t, k, x)->[0.,0.,0.], [0.])

Base.show(io::IO, cc::ControlCommand) = print(io, "(t = $(cc.t), controlfn = $(String(Symbol(cc.controlfn))), k = $(cc.k))")

# TODO: tinit not used and always passed 0
function planarcontrolfactory(prob::AttitudeProblem, tinit::Real, xinit::State, xtarget::State)
    function planarcontrol(t::Real, k::Vector, x::State)::Vector
        q = fromscalarlast(normalize(x.q))
        rotm = QuatRotation(q)
        B_B = rotm * prob.B(t)

        xI_B = rotm * [1.0, 0.0, 0.0]
        Cx_B = normalize((B_B × xI_B) × B_B)
        Cy_B = normalize(B_B × Cx_B)
        m_B = prob.m_max * (k[1] * Cx_B + k[2] * Cy_B)
        L_B = m_B × B_B
        return L_B
    end
    return planarcontrol
end

function fullcontrolfactory(prob::AttitudeProblem, tinit::Real, xinit::State, xtarget::State)
    function fullcontrol(t::Real, k::Vector, x::State)::Vector
        return k
    end
    return fullcontrol
end

function planarpdcontrolfactory(prob::AttitudeProblem, tinit::Real, xinit::State, xtarget::State)
    q_FI = fromscalarlast(normalize(xtarget.q))
    function planarpdcontrol(t::Real, k::Vector, x::State)
        q_BI = fromscalarlast(normalize(x.q))
        q_BF = q_BI*inv(q_FI)
        y = [q_BF.v1; q_BF.v2; q_BF.v3; x.ω-xtarget.ω]
        T = -[k[1]*I(3) k[2]*I(3)]*y
        rotm = QuatRotation(q_BI)
        B_B = rotm * prob.B(t)

        T = T - T⋅normalize(B_B) * normalize(B_B)
        Tmax = norm(B_B)*prob.m_max
        return norm(T) < Tmax ? T : normalize(T)*Tmax
    end
    return planarpdcontrol
end

# System Dynamics
struct ControlledAttitudeProblem
    prob::AttitudeProblem
    control::ControlCommand
end

function eulereqns(dx, x, params, t)
    """
    Derivative of Euler equations with intertially fixed dipole in constant B field
    """
    J = params.prob.J
    Jinv = params.prob.Jinv
    B_I = params.prob.B(t)
    q = fromscalarlast(normalize(x.q))

    L_B = params.control.controlfn(t, params.control.k, x)

    α = Jinv * (L_B - x.ω × (J * x.ω))
    qdot = 1 / 2 * Ξ(q) * x.ω

    dx.q = qdot
    dx.ω = α
    return nothing
end

function stepstate(state, control::ControlCommand, prob::AttitudeProblem, t0=0.0)
    params = ControlledAttitudeProblem(prob, control)
    ode = ODEProblem(eulereqns, state, t0.+(0, control.t), params)
    sol = solve(ode, reltol = 1e-7, abstol = 1e-9) # save_everystep = false
    for x in sol.u
        if norm(x.ω) > prob.ωmax
            return sol.u[end], false
        end
        if length(prob.keepouts) > 0
            rotm = QuatRotation(fromscalarlast(normalize(x.q)))
            for keepout in prob.keepouts
                obs_B = rotm * keepout.obstacle
                if angle(keepout.sensor, obs_B) < keepout.halfangle
                    return sol.u[end], false
                end
            end
        end
    end
    # @show sol.u[end]
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

timecost(roadmap, controls::Vector{ControlCommand}, istate::Integer) =
    sum([control.t for control in controls[pathto(roadmap, istate)]])

timecost(controls::Vector{ControlCommand}, path) = sum([control.t for control in controls[path]])

# Search Functions
function bonsairrt(
    prob,
    controlfnfactory::Function,
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
    kmod = nothing,
    globalbranchimprovement = false,
    timeout=Inf
)
    timeoutstart = time()
    xstart = statestart(prob)
    xgoal = stategoal(prob)

    states = [xstart]
    controls = [ControlCommand()]
    path = []
    roadmap = Dict()
    imin = 1 # Index of best state

    for iter = 1:n
        if mod(iter, 1000) == 0
            @info "on step" iter
        end
        if time() - timeoutstart > timeout
            break
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
                knew = samplecontrol()
                if !isnothing(kmod)
                    knew = controls[branch[jelbow]].k.+ kmod .* randinunit(length(controls[branch[jelbow]].k))
                    knew = norm(knew) > 1 ? knew / norm(knew) : knew # TODO generalize this
                end
                if !isnothing(tmod)
                    tnew = controls[branch[jelbow]].t * (1 - tmod + 2 * rand() * tmod)
                end
                controlfn = controlfnfactory(prob, 0, states[branch[jelbow-1]], samplestate()) # TODO: better than sample state?
                controltest = ControlCommand(tnew, controlfn, knew)

                t0 = timecost(controls, branch[1:jelbow-1])
                xfirst, valid = stepstate(
                    states[branch[jelbow-1]],
                    controltest,
                    prob,
                    t0
                )
                if valid
                    branchstates = [xfirst]
                    for control in branchcontrols
                        t0 += control.t
                        xnext, valid = stepstate(
                            branchstates[end],
                            control,
                            prob,
                            t0
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
            controlnew = nothing
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
                controlfn = controlfnfactory(prob, 0, xnear, xrand)
                t0 = timecost(roadmap, controls, inear)
                for i = 1:kbias
                    controltest = ControlCommand(sampletime(), controlfn, samplecontrol())
                    xtest, valid = stepstate(xnear, controltest, prob, t0)
                    dtest = dist(xtest, xgoal)
                    if valid & (dist(xtest, xgoal) < dcurrent)
                        xnew = xtest
                        controlnew = controltest
                        dcurrent = dtest
                    end
                end
            else
                xrand = samplestate()
                inear = argmin([dist(x, xrand) for x in states])
                xnear = states[inear]

                controlfn = controlfnfactory(prob, 0, xnear, xrand)
                controlnew = ControlCommand(sampletime(), controlfn, samplecontrol())
                t0 = timecost(roadmap, controls, inear)
                xnew, valid = stepstate(xnear, controlnew, prob, t0)
            end

            if valid
                push!(states, xnew)
                push!(controls, controlnew)
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
function lumpedcontrolfnfactory(path, controls)
    tcontrols = cumsum([control.t for control in controls[path]])
    tmax = tcontrols[end]
    function ufun(t::Real, k::Vector, x::State)
        if t == 0.0
            return controls[path[2]].controlfn(t, controls[path[2]].k, x)
        end
        icontrol = findfirst(t .< tcontrols)
        if isnothing(icontrol)
            return [0.0, 0.0, 0.0]
        end
        return controls[path[icontrol]].controlfn(t, controls[path[icontrol]].k, x)
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
    params = ControlledAttitudeProblem(prob, ControlCommand(tmax, ufun, []))
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
