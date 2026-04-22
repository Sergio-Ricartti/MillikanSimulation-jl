module Data_Process

export data_process
export plotter
export plot_vels
export regression
using JLD2
using Plots
using Statistics

include("../src/Millikan.jl")
using .Millikan

    DIRECTORY = joinpath(@__DIR__, "..", "data")

function data_process(ite::Integer)

(tray,particle,env1,dt_steps) = Read_data(ite)

end
function plot_vels(vels::Vector{Vector{Float64}}, time:: Vector{Vector{Float64}},  radius:: Vector{Float64})
    
    name = "velocity_vs_time.png"
    colores = [:crimson,:aquamarine4,:darkmagenta,:gray16]
    plt = plot(time[1], vels[1],
                title = "Velocidades de caída de tres partículas diferentes",
                ylabel = "Velocidad de caída",
                xlabel = " Tiempo (s)",
                color = colores[1],
                label = "Radio = $(radius[1])",
                legend = :bottomright)
    for i in 2:3
        plot!(time[i], vels[i],
                color = colores[i],
                label = "Radio = $(radius[i])")
    end
    scatter!([time[1][end]], [vels[1][end]],
                color = colores[4],
                marker = :circle,
                ms = 2,
                label = "Velocidad Terminal")
    for i in 2:3
        scatter!([time[i][end]], [vels[i][end]],
                    color = colores[4],
                    marker = :circle,
                    ms = 2,
                    label = "")
    end
    filepath = joinpath(DIRECTORY,"images", name)
    savefig(plt, filepath)
    name = "velocity_vs_radius.png"
    colores = [:lightslateblue, :dodgerblue2]
    vels_t = [v[end] for v in vels]
    Drho = 8.86e2-1.29
    g = 9.81
    eta = 1.81e-5
    r_teo = range(minimum(radius), maximum(radius), length=200)
    vel_teo = (2 .* (Drho) .* g .* r_teo.^2) ./ (9 .* eta)
    plt = plot(r_teo,vel_teo,
            color = colores[2],
            lw = 2,
            label = "Vf = 2r²Δρg/9η")
    scatter!([radius[1]], [vels[1][end]],
                title = "Velocidad terminal contra radio",
                ylabel = "Velocidad terminal",
                xlabel = "Radio",
                color = colores[1],
                marker = :circle,
                ms = 3,
                label = "Datos simulados N=40",
                legend = :topleft)
    for i in 2:40
        scatter!([radius[i]], [vels[i][end]],
                    color = colores[1],
                    marker = :circle,
                    ms = 3,
                    label = "")
    end
    filepath = joinpath(DIRECTORY,"images", name)
    savefig(plt,filepath)
end

function plotter(charges)
    y = 1:size(charges,1)
    name = "millikan_sim_result.png"
    
    plt = scatter(charges, y, 
                title = "Cargas de particulas", 
                marker_z = charges, 
                color = :viridis, 
                ylabel = "Número de gota", 
                xlabel = "Carga Asignada (n)",)

    filepath = joinpath(DIRECTORY,"images", name)
    savefig(plt, filepath)
end
function regression(charges::Vector{Float64})
    step = length(charges) ÷ 40  
    charges = charges[1:step:end] # Se tomarán solamente 40 cargas de 1000 (lo cual se hizo simplemente para graficar)
    charges = sort(charges) # Se ordenan los datos de menor a mayor 
    diffs = diff(charges)
    diffs_filtered = filter(d -> d >= 0.5e-19, diffs) # Se descartan diferencias menores a 0.5e-19
    e_est = median(diffs_filtered) # La mediana nos dará un valor de estimación principal para la carga elemental
    println("Valor de carga estimada: ", e_est, " C")
    n = round.(charges ./ e_est) # Se le asigna una n a cada gota

    # Regresión lineal forzada al origen (q = e*n)
    e_exp = sum(charges .* n)/sum(n .^ 2) # Solución analítica: e = sum(q*n) / sum(n^2)
    println("Valor simulado de e: ", e_exp, " C")
    println("Error %: ", round(abs(e_exp - 1.602e-19)/1.602e-19 * 100,digits = 3), " %")
    # Gráfica
    name = "Line_regression.png"
    plt = scatter(n,charges,
                    title = "Valor de e simulado por regresión lineal",
                    xlabel = "Múltiplo asignado (n)",
                    ylabel = "Carga teórica (C)",
                    label = "Cargas simuladas",
                    color = :aquamarine4,
                    marker = :circle,
                    ms = 4,
                    legend = :topleft)
    n_range = 1:maximum(n)
    e_exp = round(e_exp,sigdigits=4)
    plot!(plt, n_range, e_exp .* n_range,
            label = "e = $e_exp C",
            color = :crimson,
            linewidth = 2)
    filepath = joinpath(DIRECTORY,"images", name)
    savefig(plt, filepath)
end
function get_data(filepath::String, retrieve::String)
    data = jldopen(filepath, "r") do file
    return file[retrieve]
    end
    return data
end

function ReadData(num::Int)
    number = num
    filepath = joinpath(DIRECTORY,"millikan_test_$number.jld2")
    tray = get_data(filepath, "tray_data")
    particle = get_data(filepath, "particle_data")
    env1 = get_data(filepath, "environment_data")
    dt_steps = get_data(filepath, "time_steps")
    total_steps = size(tray,2)
    println("Success! Extracted Y-positions for ", total_steps, " steps.")
    Millikan.getInfo(particle)
    #=for (i, step_data) in enumerate(eachcol(tray))
            println(dt_steps[i], " sec | Pos: ", step_data)
    end=#
    return tray, particle, env1, dt_steps

end


end