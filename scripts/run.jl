
include("../src/Millikan.jl")
using .Millikan
include("../src/Process.jl")
using .Data_Process

using Random
using Revise
using Distributions

function run_test()
number_part = 1000
charges_acum = zeros(number_part)
tracked_part = [1,50,100]
tracked_vels = Vector{Vector{Float64}}()
tracked_time = Vector{Vector{Float64}}()
tracked_r = Vector{Float64}()
println("Running Sim...")
for n in 1:number_part
println("\n---------------------")
println("Data particle $n: ")
#GENERATE PARTICLE
#initial position
#charge generator
charge::Float64 = 1.6e-19 * rand(1:10)
#Radius distribution
rad_dis = Normal(1,0.3)
radius_rand = rand(rad_dis)*10^-6
Oil_drop = Oil_Particle(r = radius_rand, chargeQ = charge)
Env = Environment()
(trayectory, time_st, vel_hist)= runSim_Efield_Off(Oil_drop,Env, 1e-6, 1000,1e-4)
#saveData(trayectory, Oil_drop, Env, time_st)
runSim_Efield_on(Oil_drop,Env, 1e-6, 1000,1e-4, 500.00)

if n in tracked_part
    mask = time_st .> 0
    vel_y = vel_hist[2,mask] * -1
    t = time_st[mask]
    push!(tracked_vels, vel_y)
    push!(tracked_time, t)
    push!(tracked_r, radius_rand)
end

ErrorP_radius = calcErrorPorc(Oil_drop.r, Oil_drop.radius_theoric) *100
ErrorP_charge = calcErrorPorc(Oil_drop.chargeQ, Oil_drop.charge_theoric) *100

println("Terminal velocity E = 0: ", Oil_drop.velocity_t)
println("Terminal velocity E = 500V: ", Oil_drop.velocity_s)
println("Real Radius: ", Oil_drop.r, " | Exp Radius: ", Oil_drop.radius_theoric, " | % =",ErrorP_radius,"%")
println("Real Charge: ", Oil_drop.chargeQ, " | Exp Charge: ", Oil_drop.charge_theoric," | % =",ErrorP_charge,"%")
println("--------------------------------------------")
#data_process(n)
charges_acum[n] = Oil_drop.charge_theoric/Env.charge_val
end

plot_vels(tracked_vels,tracked_time,tracked_r)
plotter(charges_acum)

end

run_test()
