include("ende.jl")
include("endeplots.jl")

# SYSTEM PARAMETERS
Bmag = 0.5e-4 # Teslas (kg/As^2)
# B_hat = randunit(3)
B_hat = [0.; 0.; 1.]
B_I = B_hat * Bmag

m_max = 0.66 # A*m^2

J = [41.87 0.0 0.0; 0.0 41.87 0.0; 0.0 0.0 6.67] * 1e-3 # kg*m^2

ωmax = 5*π/180 # rad/s

prob = AttitudeProblem(B_I, m_max, J, ωmax, Quaternion(0.,0.,0.,1.), [Keepout([0.,0.,-1.],[-1.,0.,0.],50*π/180),Keepout(normalize([1.,1.,0.]),[0.,-1.,0.],20*π/180)])
# prob = AttitudeProblem(B_I, m_max, J, ωmax, Quaternion(0.,0.,0.,1.), Keepout[])
# prob = AttitudeProblem(B_I, m_max, J, ωmax)

# SOLVE

# Common Params
n = 10000
ϵ = 0.01
controlfnfactory = planarcontrolfactory
samplestate() = ComponentArray(q = toscalarlast(randquat()), ω = randinunit(3) * ωmax)
samplecontrol() = randinunit(2) # influence based on previous state?
sampletime() = rand()*25
distfn(x1, x2) = distB(x1, x2, prob.B, 30, 0, 0)

pbias = 0.66
kbias = 25
biasfactor = Inf

pbend = 0.03
kbend = 50
bendfactor = Inf
tmod = 0.1
kmod = nothing
globalbranchimprovement = true

alg_mode = "bonsai"
control_mode = "pd"
timeout = Inf

if alg_mode == "rrt"
    pbend = 0
elseif alg_mode == "bonsai"
    pbend = 0.03
end

if control_mode == "const"
    controlfnfactory = planarcontrolfactory
    samplecontrol() = randinunit(2)
elseif control_mode == "pd"
    controlfnfactory = planarpdcontrolfactory
    samplecontrol() = [rand()*0.005; rand()*0.03]
end


path, roadmap, states, controls = bonsairrt(prob, controlfnfactory, samplecontrol, samplestate, sampletime, distfn, n, ϵ, pbias, kbias, biasfactor, pbend, kbend, bendfactor;
                                            tmod=tmod, kmod=kmod, globalbranchimprovement=globalbranchimprovement, timeout=timeout)

suffix = string(alg_mode,"-",control_mode)

# POSTPROCESSING
framerate = 10
timefactor = 10
saveat = timefactor/framerate
comparison = states[path[end]]

(ufun, tmax, tstops) = lumpedcontrolfnfactory(path, controls)
sol = resamplesolution(prob, ufun, tmax, saveat=saveat, comparison=comparison, tstops=tstops)

# PLOTTING

CairoMakie.activate!()

fig = Figure(resolution=(2.25,3).*72.0.*2.5, font = "CMU Serif", fontsize = 18)
ax1 = Axis(fig[1,1])
roadmapplot2d(ax1, path, roadmap, states, controls, prob, distfn)
# ax2 = Axis3(fig[1,2])
# roadmapplot3d(ax2, path, roadmap, states, controls, prob)
# ax2.azimuth = -1.25
# ax2.elevation = 0.5
save("output/roadmap_$(suffix).pdf", fig, pt_per_unit=1/2.5)


# saveas = "output/time_animation_$(suffix).mp4"
# saveanimation(sol, ufun, prob, saveas, framerate)


fig = Figure(resolution=(4,3).*72.0.*2.5, font = "CMU Serif", fontsize = 18)
ax1 = Axis(fig[2,1])
[lines!(ax1, sol.t, [(QuatRotation(normalize(inv(fromscalarlast(x.q))))*x.ω)[i]*180/π for x in sol.u], label=l) for (i, l) in enumerate([L"$\omega_x$", L"$\omega_y$", L"$\omega_z$"])]
ax1.ylabel = L"$ω$ [deg/s]"
fig[2,2] = Legend(fig, ax1)

ax2 = Axis(fig[1,1])
roll = [atand(2(x.q[4]*x.q[1]+x.q[2]*x.q[3]), 1-2(x.q[1]^2+x.q[2]^2)) for x in sol.u]
pitch = [asind(2(x.q[4]*x.q[2]-x.q[3]*x.q[1])) for x in sol.u]
yaw = [atand(2(x.q[4]*x.q[3]+x.q[1]*x.q[2]), 1-2*(x.q[2]^2+x.q[3]^2)) for x in sol.u]

lines!(ax2, sol.t, roll, label=L"$\phi$")
lines!(ax2, sol.t, pitch, label=L"$\theta$")
lines!(ax2, sol.t, yaw, label=L"$\psi$")
lines!(ax2, sol.t, [2*acosd(x.q[4]) for x in sol.u], label="Eig. Ang.", color=:black)
ax2.ylabel = L"Angle [deg]$ $"
ax2.yticks = -180:90:180
fig[1,2] = Legend(fig, ax2)

ax3 = Axis(fig[3,1])
# ax3.ytickformat =
[lines!(ax3, sol.t, [(QuatRotation(normalize(inv(fromscalarlast(x.q))))*ufun(t, [], x))[i] for (t,x) in zip(sol.t, sol.u)], label=l) for (i, l) in enumerate([L"$L_x$", L"$L_y$", L"$L_z$"])]
ax3.ylabel = L"$\vec{L}$ [N⋅m]"
fig[3,2] = Legend(fig, ax3)

ax3.xlabel = "Time [s]"
linkxaxes!(ax1,ax2,ax3)
save("output/solutiondynamics_$(suffix).pdf", fig, pt_per_unit=1/2.5)


fig = Figure(resolution=(2.75,2.5).*72.0.*2.5, font = "CMU Serif", fontsize = 18)
ax1 = Axis(fig[1,1], yscale=log10)
lines!(ax1, sol.t, [2*acosd(abs(x.q[4])) for x in sol.u], label="Eig. Angle")
ax1.ylabel = L"Eig. Angle [deg] $$"

ax2 = Axis(fig[2,1], yscale=log10)#, backgroundcolor = :transparent, yaxisposition = :right)
lines!(ax2, sol.t[2:end], [norm(x.ω)*180/π for x in sol.u[2:end]])
ax2.ylabel = L"$|\vec{\omega}|$ [deg/s]"

ax2.xlabel = "Time [s]"
linkxaxes!(ax1,ax2)
save("output/convergence_$(suffix).pdf", fig, pt_per_unit=1/2.5)


fig = Figure(resolution=(3.5,2.5).*72.0.*2.5, font = "CMU Serif", fontsize = 18)
ax = Axis(fig[1,1])
latlongkeepout(ax, sol, prob)
save("output/latlong_$(suffix).pdf", fig, pt_per_unit=1/2.5)
