# load in the data

module loadData

export load_plot_data

using JLD2
using Oceananigans

function load_plot_data(name, numSizeClasses)

    n_dict = Dict(Symbol("n$n") => FieldTimeSeries(name, "n$n") for n in 1:numSizeClasses)
    # Merge with fixed fields and convert to NamedTuple
    time_series = (; 
        w = FieldTimeSeries(name, "w"),
        u = FieldTimeSeries(name, "u"),
        v = FieldTimeSeries(name, "v"),
        T = FieldTimeSeries(name, "T"),
        S = FieldTimeSeries(name, "S"),
        n_dict...)   # Expand dictionary into NamedTuple fields

    return time_series
end

end
