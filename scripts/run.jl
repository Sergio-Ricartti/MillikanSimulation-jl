
include("../src/Millikan.jl")
using .Millikan
include("../src/Process.jl")
using .Data_Process

using Random
using Revise
using Distributions

function run_test()
number_part = 10
charges_acum = zeros(number_part)
println("Running Sim...")
for n in 1:number_part
println("\n---------------------")
println("Data particle $n: ")
#GENERATE PARTICLE
#initial position
ini_pos = [0.0,1e-5,0.0]
#charge generator
charge::Float64 = 1.6e-19 * rand(1:5)
#Radius distribution
rad_dis = Normal(1,0.3)
radius_rand = rand(rad_dis)*10^-6
Oil_drop = Oil_Particle(r = radius_rand, chargeQ = charge, pos = ini_pos)
Env = Environment()
(trayectory, time_st)= runSim(Oil_drop,Env, 1e-6, 1000,1e-3)
saveData(trayectory, Oil_drop, Env, time_st)
ErrorP_radius = calcErrorPorc(Oil_drop.r, Oil_drop.radius_theoric) *100

println("Terminal velocity E = 0: ", Oil_drop.velocity_t)
println("Real Radius: ", Oil_drop.r, " | Exp Radius: ", Oil_drop.radius_theoric, " | % =",ErrorP_radius,"%")
println("Real Charge: ", Oil_drop.chargeQ, " | Exp Charge: ", Oil_drop.charge_theoric)
#data_process(n)

end



end

run_test()
