module Millikan
export Oil_Particle
export Environment
export runSim_Efield_Off
export runSim_Efield_on

export getInfo
export saveData
export calcErrorPorc

using LinearAlgebra
using JLD2



Base.@kwdef mutable struct Oil_Particle        #Caracteristicas de la particula
r::Float64                       #m^2
mass::Float64 = 0.0              #kg    - Depends on radius r
pos::Vector{Float64} = [0.0,1e-5,0.0]#m
vel::Vector{Float64} = [0,0,0]   #m/s
acc::Vector{Float64} = [0,0,0]   #m/s^2
chargeQ::Float64                 #C
volume::Float64 = 0.0            #m^3
velocity_t::Vector{Float64} = [0,0,0]   #m/s  - Terminal velocity (E=0)
velocity_s::Vector{Float64} = [0,0,0]   #m/s  -     
charge_theoric::Float64 = 0.0
radius_theoric::Float64 = 0.0
voltage::Float64 = 0.0
end

function getInfo(p::Oil_Particle)
println("||")
println("Mass:", p.mass, "kg")
println("Volume:", p.volume, "m^3")
println("Charge:", p.chargeQ, "C")
println("Radius:", p.r, "m")
println("||")
end

# Valores de variables constantes o del entorno
Base.@kwdef struct Environment
g::Vector{Float64} = [0,-9.81,0]            #m/s^2
eta::Float64 = 1.81e-5                      #Pa*s
rho_oil::Float64 = 8.86e2                   #kg/m^3
rho_air::Float64 = 1.29                     #kg/m^3
charge_val::Float64 = 1.6e-19               #C
distance_d::Vector{Float64} = [0,1e-4,0]    #m - Distance between plates is 0.1mm (Millikan., 1913)
end

#metodo iterativo de euler
function Euler!(p::Oil_Particle,env::Environment, dt::Float64)
    println("-----------------------------")
    Forces = ForcesCalc(p,env, p.vel)
    p.acc .= Forces / p.mass
    println("Acceleration = ", p.acc)
    p.vel .+= p.acc .* dt
    p.pos .+= p.vel .* dt
end
#Metodo iterativo RK4
function RK4!(p::Oil_Particle,env::Environment, dt::Float64)
#Stage1
    vel1 = p.vel
    acc1 = ForcesCalc(p,env, vel1)/p.mass
#Stage 2
    vel2 = p.vel .+dt .* acc1/2
    acc2 = ForcesCalc(p,env,vel2)/p.mass
#Stage 3
    vel3 = p.vel .+dt .*acc2/2
    acc3 = ForcesCalc(p,env,vel3)/p.mass
#Stage 4
    vel4 = p.vel .+ dt .*acc3
    acc4 = ForcesCalc(p, env, vel4)/p.mass
#Final Values
    p.pos .+= (dt/6) .* (vel1 .+ 2 .* vel2 .+ 2 .* vel3 .+ vel4)
    p.vel .+= (dt / 6) .* (acc1 .+ 2 .* acc2 .+ 2 .* acc3 .+ acc4)
    p.acc .= ForcesCalc(p,env,p.vel) /p.mass

    #println("Acceleration: ", p.acc)
    #println("Velocity: ", p.vel)
end

#Calculo de fuerzas
#Fuerza de gravedad effectiva
function GForce(p::Oil_Particle, env::Environment)
    return env.g*(p.mass - p.volume * env.rho_air)
end

#Fuerza de arrastre
function DragForce(p::Oil_Particle, env::Environment, velocity::Vector)
    return -6*env.eta*pi*p.r *velocity
end
#Fuerza Electrica
function ElectricForce(p::Oil_Particle, env::Environment)
    d_mag = norm(env.distance_d)
    F_e = p.chargeQ .*(p.voltage/d_mag) .* (env.distance_d ./d_mag)
    return F_e
end
function ForcesCalc(p::Oil_Particle, env::Environment, vel::Vector)
    Fg = GForce(p,env)
    Fd = DragForce(p,env, vel)
    Fe = ElectricForce(p,env)
    F_total = Fg .+ Fd .+ Fe
    return F_total

end


function saveData(traj::Array, p::Oil_Particle, env::Environment, time::Array)

DIRECTORY = joinpath(@__DIR__, "..", "data")
name = "millikan_test_"
i = 1
while isfile(joinpath(DIRECTORY,"$(name)$(i).jld2"))
    i+=1
end
full_name = "$name$i"
full_path = joinpath(DIRECTORY,"$full_name.jld2")
jldsave(full_path;
    tray_data = traj,
    particle_data = p,
    environment_data = env,
    time_steps = time)
#println("Data succesfully saved to -> ", full_path)

end

function calcErrorPorc(t, e)
    return abs((t-e)/t)
end



function runSim_Efield_Off(p::Oil_Particle,env::Environment, dt::Float64, steps::Int, error::Float64)
    p.volume = (4/3)*pi * p.r^3
    p.mass = p.volume * env.rho_oil
    vel_error = 1.0
    p.voltage = 0 #V
    p.pos = [0.0,1e-5,0.0]
    p.vel = [0,0,0]
    p.acc = [0,0,0]
    getInfo(p)
    pos_history = zeros(Float64, 3, steps)
    vel_history = zeros(Float64, 3, steps)
    time_steps = zeros(Float64,steps)
    dt_accum = 0.0
    i = 1
    while i<=steps && vel_error>=error      #Obtener velocidad terminal realizando la simulacion con V = 0
        past_vel_y = p.vel[2]
        RK4!(p,env, dt)
        vel_error = calcErrorPorc(past_vel_y, p.vel[2])
        pos_history[:,i] .= p.pos
        vel_history[:,i] .= p.vel
        time_steps[i] = dt_accum
        dt_accum += dt
        #println(dt_accum, "|", p.pos)
        i +=1
    end
    println("Vel_Error> ", vel_error )

    p.velocity_t = p.vel
    p.radius_theoric = sqrt((9*env.eta*norm(p.velocity_t))/(2*norm(env.g)*(env.rho_oil - env.rho_air)))

    return pos_history, time_steps, vel_history

end

function runSim_Efield_on(p::Oil_Particle, env::Environment, dt::Float64, steps::Int, error::Float64, volt::Float64)
    vel_error = 1.0
    i = 1
    p.voltage = volt #V
    p.pos = [0.0,1e-5,0.0]
    p.vel = [0,0,0]
    p.acc = [0,0,0]

    while i<=steps && vel_error>=error      #Obtener velocidad terminal realizando la simulacion con V = 0
            past_vel_y = p.vel[2]
            RK4!(p,env, dt)
            vel_error = calcErrorPorc(past_vel_y, p.vel[2])
            #println(dt_accum, "|", p.pos)
            i +=1
    end
    println("Vel_Error> ", vel_error )
    p.velocity_s = p.vel
    p.charge_theoric = (6*pi*p.radius_theoric*env.eta*norm(env.distance_d)/(p.voltage))*(norm(p.velocity_t)+norm(p.velocity_s))

end



end