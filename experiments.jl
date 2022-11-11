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
n = 20000
ϵ = 0.03
samplestate() = ComponentArray(q = toscalarlast(randquat()), ω = randinunit(3) * ωmax)
samplecontrol() = randinunit(2) # influence based on previous state?
sampletime() = rand()*25
distfn(x1, x2) = distB(x1, x2, prob.B, 30, 0, 0)

pbias = 0.66
kbias = 25
biasfactor = 20

pbend = 0.03
kbend = 50
bendfactor = 20
tmod = 0.1
umod = nothing
globalbranchimprovement = true

standardrrt = true
if standardrrt
    pbend = 0.
    # biasfactor = Inf
end
suffix = "biasrrt"

path, roadmap, states, controls = bonsairrt(prob, samplecontrol, samplestate, sampletime, distfn, n, ϵ, pbias, kbias, biasfactor, pbend, kbend, bendfactor; tmod=tmod, umod=umod, globalbranchimprovement=globalbranchimprovement)


# δₛ = 0.05
# δbn = 2δₛ
# path, roadmap, states, controls = sst(prob, samplecontrol, samplestate, sampletime, distfn, n, ϵ, δbn, δₛ, pbias)

# POSTPROCESSING
framerate = 10
timefactor = 10
saveat = timefactor/framerate
comparison = states[path[end]]

(ufun, tmax, tstops) = controlfnfactory(path, controls)
sol = resamplesolution(prob, ufun, tmax, saveat=saveat, comparison=comparison, tstops=tstops)

# PLOTTING
fig = Figure(resolution=(1280,550).*2, fontsize=24)
ax1 = Axis(fig[1,1])
roadmapplot2d(ax1, path, roadmap, states, controls, prob, distfn)
ax2 = Axis3(fig[1,2])
roadmapplot3d(ax2, path, roadmap, states, controls, prob)
ax2.azimuth = -1.25
ax2.elevation = 0.5
save("output/roadmap_$(suffix).png", fig)

saveas = "output/time_animation_$(suffix).mp4"
saveanimation(sol, ufun, prob, saveas, framerate)

fig = Figure(resolution=(800,600).*2, fontsize=24)
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
save("output/solutiondynamics_$(suffix).png", fig)
