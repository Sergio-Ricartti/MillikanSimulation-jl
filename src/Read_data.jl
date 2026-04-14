module data_READ
include("../src/Millikan.jl")

export ReadData
export plotter
using JLD2
using .Millikan
using Plots

function get_data(filepath::String, retrieve::String)
    data = jldopen(filepath, "r") do file
    return file[retrieve]
    end
    return data
end

function ReadData(num::Int)
    DIRECTORY = joinpath(@__DIR__, "..", "data")
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