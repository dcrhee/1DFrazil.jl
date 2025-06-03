# load in the data

module loadLastSnapshot

export load_data

using JLD2
using Oceananigans
using Statistics

function find_z_mean(var)
    # finds the z-mean of the array
    return vec(mean(var[1, 1, :, :], dims = 1))
end

function load_data(name, numSizeClasses)

    n_dict = Dict(Symbol("n$n") => FieldTimeSeries(name, "n$n") for n in 1:numSizeClasses)
    # Merge with fixed fields and convert to NamedTuple
    time_series = (; 
        w = FieldTimeSeries(name, "w"),
        u = FieldTimeSeries(name, "u"),
        v = FieldTimeSeries(name, "v"),
        T = FieldTimeSeries(name, "T"),
        S = FieldTimeSeries(name, "S"),
        n_dict...)   # Expand dictionary into NamedTuple fields

    end_values = zeros(numSizeClasses)
    for i in 1:numSizeClasses
        nvals = find_z_mean(getfield(time_series, Symbol("n$i")))
        end_values[i] = nvals[end]
    end

    return end_values
end

end
