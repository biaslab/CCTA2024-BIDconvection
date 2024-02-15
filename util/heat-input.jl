using LinearAlgebra


function u(t::Float64, b::Float64; i::Int64=1) 
    "Constant heater"
    return b
end

function u(t::Float64, b::Float64, t_on::Float64; i::Int64=1, β::Float64=0.1) 
    "Turn on heater"
    return b ./(1+exp(-β*(t - t_on))) 
end

function u(t::Float64; i::Int64=1, b::Float64=1.0, t_on::Float64=20., t_off::Float64=60., β1::Float64=0.6, β2::Float64=0.4) 
    "Pulse heater"
    return b ./(1+exp(-β1*(t - t_on))) * exp(-β2*(t - t_off))./(1+exp(-β2*(t - t_off)))
end