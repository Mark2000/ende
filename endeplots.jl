using Makie.GeometryBasics
using GLMakie
using FileIO

Axis = GLMakie.Axis

include("ende.jl")

function roadmapplot2d(ax::Axis, path, roadmap, states, controls, prob::AttitudeProblem, distfn::Function)
    t_rm = [timecost(roadmap, controls, p) for pair in roadmap for p in pair]
    d_rm = [distfn(states[p], stategoal(prob)) for pair in roadmap for p in pair]

    t_path = [timecost(roadmap, controls, p) for p in path]
    d_path = [distfn(states[p], stategoal(prob)) for p in path]

    linesegments!(ax, d_rm, t_rm, color=:lightgrey)
    lines!(ax, d_path, t_path, color=:red)
    scatter!(ax, [d_path[1]], [0], color=:red, markersize=20)
    lines!(ax, [0,0], [0, max(t_rm...)], color=:green)
    ax.xlabel = "Dist. to Goal"
    ax.ylabel = "Path Length [s]"
    return ax
end

function roadmapplot3d(ax::Axis3, path, roadmap, states, controls, prob::AttitudeProblem)
    t_rm = [timecost(roadmap, controls, p) for pair in roadmap for p in pair]
    ω_rm = [dist(states[p], stategoal(prob), 1e10)/1e10 for pair in roadmap for p in pair]
    ang_rm = [dist(states[p], stategoal(prob), 0) for pair in roadmap for p in pair]

    t_path = [timecost(roadmap, controls, p) for p in path]
    ω_path = [dist(states[p], stategoal(prob), 1e10)/1e10 for p in path]
    ang_path = [dist(states[p], stategoal(prob), 0) for p in path]

    linesegments!(ax, t_rm, ω_rm, ang_rm, color=(:lightgrey, 0.5), transparency=true)
    lines!(ax, t_path, ω_path, ang_path, color=:red)
    lines!(ax, [0, max(t_rm...)], [0,0], [0,0], color=:green)
    ax.xlabel = "Path Length [s]"
    ax.ylabel = "Angle Rate Error [rad/s]"
    ax.zlabel = "Angle Error [rad]"
    return ax
end

function saveanimation(sol, ufun::Function, prob::AttitudeProblem, path, framerate, resolution=(1280,1280))
    fig = Figure(resolution=resolution)
    sat = load("dependencies/2RU-GenericCubesat.stl")

    tstamp = Observable(1)

    ax = Axis3(fig[1,1], aspect = :data, title = @lift("t = $(round(sol.t[$tstamp], digits = 1))"))
    xlims!(ax,(-3.,3))
    ylims!(ax,(-3.,3))
    zlims!(ax,(-3.,3))

    B_hat = normalized(prob.B)

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

    record(fig, path, timestamps;
            framerate = framerate) do t
        tstamp[] = t
    end
    @info "Animation saved!"
end
