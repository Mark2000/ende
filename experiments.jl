include("ende.jl")
include("endeplots.jl")

# SYSTEM PARAMETERS
Bmag = 0.5e-4 # Teslas (kg/As^2)
# B_hat = randunit(3)
B_hat = [0.; 0.; 1.]
B_I = B_hat * Bmag

m_max = 0.66 # A*m^2

J = [41.87 0.0 0.0; 0.0 41.87 0.0; 0.0 0.0 6.67] * 1e-3 # kg*m^2

ωmax = 2*π/180 # rad/s

prob = AttitudeProblem(B_I, m_max, J, ωmax, Quaternion(0.,0.,0.,1.), [Keepout([0.,0.,-1.],[-1.,0.,0.],50*π/180),Keepout(normalize([1.,1.,0.]),[0.,-1.,0.],20*π/180)])
# prob = AttitudeProblem(B_I, m_max, J, ωmax, Quaternion(0.,0.,0.,1.), Keepout[])
# prob = AttitudeProblem(B_I, m_max, J, ωmax)

# SOLVE

# Common Params
n = 15000
ϵ = 0.10
pbias = 0.10
samplestate() = ComponentArray(q = toscalarlast(randquat()), ω = randinunit(3) * ωmax)
samplecontrol() = randinunit(2)
sampletime() = rand()*10
distfn(x1, x2) = distB(x1, x2, prob.B, 20, 0, 0)

k = 25
path, roadmap, states, controls = dynamicrrt(prob, samplecontrol, samplestate, sampletime, distfn, n, ϵ, pbias, k)

# δₛ = 0.025
# δbn = 2δₛ
# path, roadmap, states, controls = sst(prob, samplecontrol, samplestate, sampletime, distfn, n, ϵ, δbn, δₛ, pbias)

# POSTPROCESSING
framerate = 10
timefactor = 10
saveat = timefactor/framerate
comparison = states[path[end]]

(ufun, tmax) = controlfnfactory(path, controls)
sol = resamplesolution(prob, ufun, tmax, saveat=saveat, comparison=comparison)

# Plotting
path = "results/time_animation.mp4"
saveanimation(sol, ufun, prob, path, framerate)

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
save("results/solutiondynamics.png", fig)

fig = Figure()
ax1 = Axis(fig[1,1])
roadmapplot2d(ax1, path, roadmap, states, controls, prob, distfn)
ax2 = Axis3(fig[1,2])
roadmapplot3d(ax2, path, roadmap, states, controls, prob)
save("results/roadmaps.png", fig)
