module Data_Process
export data_process
include("../src/Read_data.jl")
using .data_READ
include("../src/Millikan.jl")
using .Millikan

function data_process(ite::Integer)

(tray,particle,env1,dt_steps) = Read_data(ite)


end

function plotter(trajectory::Array, timesteps::Array)
    a = 1
    
end



end