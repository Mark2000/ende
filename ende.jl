using ComponentArrays
using DifferentialEquations
using LinearAlgebra
using Optim
using Quaternions
using Random
using Rotations

using Makie.GeometryBasics
using GLMakie
using FileIO

Quaternion = Quaternions.Quaternion
Axis = GLMakie.Axis

randquat() = normalize(Quaternion(randn(4)...))

randunit(n = 3) = normalize([randn(n)...])

function randinunit(n = 2)
    while true
        v = 2*rand(n).-1
        if norm(v) ≤ 1
            return v
        end
    end
end

fromscalarlast(q) = Quaternion(q[4], q[1], q[2], q[3])
toscalarlast(q) = [q.v1; q.v2; q.v3; q.s]

Ξ(q) = [q.s -q.v3 q.v2; q.v3 q.s -q.v1; -q.v2 q.v1 q.s; -q.v1 -q.v2 -q.v3]

function dist(x1, x2, α = 5)
    # dq = conj(fromscalarlast(x1.q)) * fromscalarlast(x2.q)
    # dangle = abs(rotation_angle(QuatRotation(normalize(dq))))
    qw = max(min(x1.q⋅x2.q,1.),-1.)
    dangle = abs(2*acos(qw))

    dω = norm(x1.ω - x2.ω)

    dangle + α * dω
end

function distB(x1, x2, B, α=5, ϵ=1, β=1)
    # dq = conj(fromscalarlast(x1.q)) * fromscalarlast(x2.q)
    # dangle = abs(rotation_angle(QuatRotation(normalize(dq))))
    if ϵ==β==0
        return dist(x1, x2, α)
    end

    B_hat = normalize(B)
    qw = max(min(x1.q⋅x2.q,1.),-1.)
    dangle = abs(2*acos(qw)) # in range [0,π]
    qv = -x1.q[1:3].*x2.q[4] + x1.q[4].*x2.q[1:3] - x1.q[1:3]×x2.q[1:3]
    eigax = qv/sqrt(1-qw^2)
    oopangle = abs(B_hat⋅eigax) # in range [0,1]

    dω = norm(x1.ω - x2.ω) # in range [0,ωmax]
    oopω = abs(B_hat⋅normalize(x1.ω-x2.ω))  # in range [0,1]
    if isnan(oopω)
        oopω = 0
    end

    cost = dangle * (1 + ϵ*oopangle) + α * dω * (1 + β*oopω)

    # @info cost dangle oopangle dω oopω
    return cost
end

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

# TODO include Bhat
AttitudeProblem(B::Vector, m_max::Real, J::Matrix, ωmax::Real) = AttitudeProblem(B, m_max, J, inv(J), randquat(), [0.0; 0.0; 0.0], Quaternion(1.), [0.0; 0.0; 0.0], ωmax, Keepout[])
AttitudeProblem(B::Vector, m_max::Real, J::Matrix, ωmax::Real, qstart::Quaternion) = AttitudeProblem(B, m_max, J, inv(J), qstart, [0.0; 0.0; 0.0], Quaternion(1.), [0.0; 0.0; 0.0], ωmax, Keepout[])
AttitudeProblem(B::Vector, m_max::Real, J::Matrix, ωmax::Real, qstart::Quaternion, keepouts::Vector{Keepout}) = AttitudeProblem(B, m_max, J, inv(J), qstart, [0.0; 0.0; 0.0], Quaternion(1.), [0.0; 0.0; 0.0], ωmax, keepouts)

statestart(prob::AttitudeProblem) = ComponentArray(q = toscalarlast(prob.qstart), ω = prob.ωstart)
stategoal(prob::AttitudeProblem) = ComponentArray(q = toscalarlast(prob.qgoal), ω = prob.ωgoal)


struct ControlledAttitudeProblem
    prob::AttitudeProblem
    u::Function
end

function eulereqns(dx, x, params, t)
    """
    Derivative of Euler equations with intertially fixed dipole in constant B field
    """
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

function angle(a, b)
    return acos(clamp(a⋅b/(norm(a)*norm(b)), -1, 1))
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
                obs_B = rotm*keepout.obstacle
                if angle(keepout.sensor, obs_B) < keepout.halfangle
                    return sol.u[end], false
                end
            end
        end
    end
    sol.u[end], true
end

function pathto(roadmap, istate)
    path = [istate]
    while path[1] != 1
        pushfirst!(path, roadmap[path[1]])
    end
    return path
end

struct ControlCommand
    t::Real
    u::Vector
end

function dynamicrrt(prob, samplecontrol::Function, samplestate::Function, sampletime::Function, dist::Function, n, ϵ, pbias, k)
    xstart = statestart(prob)
    xgoal = stategoal(prob)

    states = [xstart]
    controls = [ControlCommand(0., [0.0; 0.0])]
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
        dt = 0.

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

# SST FUNCTIONS
function sst(prob, samplecontrol::Function, samplestate::Function, sampletime::Function, dist::Function, n, ϵ, δbn, δₛ, pbias)
    # Currently assumes time is cost
    xstart = statestart(prob)
    xgoal = stategoal(prob)

    states = [xstart]
    controls = [ControlCommand(0., [0.0; 0.0])]
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
        isel = isempty(Inear) ? Vactive[argmin([dist(x, xrand) for x in states[Vactive]])] : Inear[argmin([cost(i) for i in Inear])]
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
            isnew = collect(keys(S))[argmin([dist(x, xnew) for x in states[collect(keys(S))]])]
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

                if dist(xnew, xgoal) < dist(states[imin], xgoal) || ((dist(xnew, xgoal) < ϵ) && (cost(ixnew) < cost(imin)))
                    imin = length(states)
                    if xgoal == xrand
                        @info "Biased step"
                    end
                    @info xnew imin dist(xnew, xgoal) cost(imin)
                end

                # Prune dominated nodes
                if !isnothing(ixpeer)
                    filter!(x->x≠ixpeer,Vactive)
                    push!(Vinactive, ixpeer)
                end
                S[isnew] = ixnew
                while !isnothing(ixpeer) && ixpeer ∉ values(roadmap) && ixpeer ∈ Vinactive
                    ixparent = roadmap[ixpeer]
                    delete!(roadmap, ixpeer)
                    filter!(x->x≠ixpeer,Vinactive)
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

# Execute Problem
Bmag = 0.5e-4 # Teslas (kg/As^2)
# B_hat = randunit()
B_hat = [0.; 0.; 1.]
B_I = B_hat * Bmag

m_max = 0.66 # A*m^2

J = [41.87 0.0 0.0; 0.0 41.87 0.0; 0.0 0.0 6.67] * 1e-3 # kg*m^2

ωmax = 20*π/180

prob = AttitudeProblem(B_I, m_max, J, ωmax, Quaternion(0.,0.,0.,1.), [Keepout([0.,0.,-1.],[-1.,0.,0.],50*π/180),Keepout(normalize([1.,1.,0.]),[0.,-1.,0.],20*π/180)])
prob = AttitudeProblem(B_I, m_max, J, ωmax, Quaternion(0.,0.,0.,1.), Keepout[])
# prob = AttitudeProblem(B_I, m_max, J, ωmax)


samplestate() = ComponentArray(q = toscalarlast(randquat()), ω = randinunit(3) * ωmax)
samplecontrol() = randinunit(2)
sampletime() = rand()*10
distfn(x1, x2) = distB(x1, x2, prob.B, 5, 0, 0)

n = 10000
ϵ = 0.10

pbias = 0.10
k = 25

# path, roadmap, states, controls = dynamicrrt(prob, samplecontrol, samplestate, sampletime, distfn, n, ϵ, pbias, k)


δₛ = 0.025
δbn = 2δₛ
path, roadmap, states, controls = sst(prob, samplecontrol, samplestate, sampletime, distfn, n, ϵ, δbn, δₛ, pbias)


# High-res solution
framerate = 10
timefactor = 10

tcontrols = cumsum([control.t for control in controls[path]])
tmax = tcontrols[end]
function ufun(t)
    if t == 0.0
        return controls[path[2]].u
    end
    icontrol = findfirst(t .< tcontrols)
    if isnothing(icontrol)
        return [0.,0.]
    end
    return controls[path[icontrol]].u
end

params = ControlledAttitudeProblem(prob, ufun)
ode = ODEProblem(eulereqns, statestart(prob), (0, tmax), params)
sol = solve(ode, reltol = 1e-7, abstol = 1e-9, saveat=timefactor/framerate, tstops=tcontrols)
@info "Integrated difference should be small" dist(sol.u[end], states[path[end]])

# Plotting
tstamp = Observable(1)
sat = load("2RU-GenericCubesat.stl")

fig = Figure(resolution=(1280,1280))
ax = Axis3(fig[1,1], aspect = :data, title = @lift("t = $(round(sol.t[$tstamp], digits = 1))"))
xlims!(ax,(-3.,3))
ylims!(ax,(-3.,3))
zlims!(ax,(-3.,3))

# q = normalize(fromscalarlast(xgoal.q))
# wireframe!(ax, GeometryBasics.Mesh([QuatRotation(q)*p for p in sat.position],faces(sat)), color=(:black,0.025))

arrows!(ax,[Point3([0.; 0.; 0.])], [Point3(B_hat)], color=:green, lengthscale=2, label="B")

ωvec = @lift( [Point3( (QuatRotation(normalize(inv(fromscalarlast(sol.u[$tstamp].q)))) * sol.u[$tstamp].ω)... )] )
arrows!(ax, [Point3([0.; 0.; 0.])], ωvec, color=:purple, lengthscale=20, label="ω")

Bx = normalize((B_hat × [1.0, 0.0, 0.0]) × B_hat)
By = B_hat × Bx
uvec = @lift( [Point3( ufun(sol.t[$tstamp])[1].*Bx + ufun(sol.t[$tstamp])[2].*By )] )
arrows!(ax, [Point3([0.; 0.; 0.])], uvec, color=:red, lengthscale=2, label="m")
Lvec = @lift( [Point3( (ufun(sol.t[$tstamp])[1].*Bx + ufun(sol.t[$tstamp])[2].*By) × B_hat)] )
ωvecscale = @lift( $ωvec.*20 )
arrows!(ax, ωvecscale, Lvec, color=:blue, lengthscale=2, label="L")

rotm = @lift( QuatRotation(normalize(inv(fromscalarlast(sol.u[$tstamp].q)))) )

keepoutcolors = [:yellow, :orange, :pink]
for (keepout, color) in zip(prob.keepouts, keepoutcolors)
    sensorvec = @lift( [$rotm*keepout.sensor] )
    arrows!(ax, [Point3([0.; 0.; 0.])], sensorvec, color=color, lengthscale=2, label="sensor")

    for θ in 0:0.03:2π
        r = 3.
        obstacle = keepout.obstacle
        halfangle = keepout.halfangle
        d = r/tan(halfangle)

        rx = abs.(obstacle) == [0.0, 0.0, 1.0] ? [1.0, 0.0, 0.0] : obstacle×[0.0, 0.0, 1.0]

        ry = obstacle×rx

        point = obstacle*d + r*cos(θ)*rx + r*sin(θ)*ry
        # @show point
        lines!(ax, [Point3([0.; 0.; 0.]), Point3(point)], color=(color,0.5))
    end
end

satrot = @lift(GeometryBasics.Mesh([$rotm*p for p in sat.position],faces(sat)))
mesh!(ax, satrot, color=(:grey, 0.2))


timestamps = 1:length(sol.t)

record(fig, "time_animation.mp4", timestamps;
        framerate = framerate) do t
    tstamp[] = t
end
@info "Animation saved!"

fig = Figure()
ax1 = Axis(fig[1,1])
[lines!(ax1, sol.t, [(QuatRotation(normalize(inv(fromscalarlast(x.q))))*x.ω)[i]*180/π for x in sol.u]) for i in 1:3]
ax1.ylabel = "ω [deg/s]"

ax2 = Axis(fig[2,1])
[lines!(ax2, sol.t, [ufun(t)[i] for t in sol.t]) for i in 1:2]
ax2.ylabel = "u [-]"

ax3 = Axis(fig[3,1])
[lines!(ax3, sol.t, [x.q[i] for x in sol.u]) for i in 1:4]
ax3.ylabel = "q [-]"

ax3.xlabel = "Time [s]"

linkxaxes!(ax1,ax2,ax3)

save("results.png", fig)



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
