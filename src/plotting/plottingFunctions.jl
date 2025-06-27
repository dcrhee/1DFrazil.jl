module plottingFunctions

export find_z_mean, find_C_const, generatePDF, generateCPDF

using Statistics
using Oceananigans
using KernelDensity

function find_z_mean(var)
    # finds the z-mean of the array
    return vec(mean(var[1, 1, :, :], dims = 1))
end

function find_C_const(indx, Rs, aspect_ratio, volume)
    # finds the z-mean of the array
    return 2*π*Rs[indx]^3/(aspect_ratio*volume)
end

# generate pdf
function generatePDF(tindx, chosentimeseries, Rs, numSizeClasses)
    # input: tindx = time
    expanded_data = Vector{Float64}()  # Pre-allocate an empty vector
    barvals = zeros(length(Rs)-1)
    for i in 2:numSizeClasses
        
        nivals = find_z_mean(getfield(chosentimeseries, Symbol("n$i")))  # Extract tracer values
        println(i, nivals[tindx])
        append!(expanded_data, Rs[i] * nivals[tindx])#ones(round(Int, nivals[tindx])))
        barvals[i-1] = nivals[tindx]
    end
        
    # Perform Kernel Density Estimation of everything apart from size class 1
    #pdfsol = kde(expanded_data)
    #pdfsol = kde((Rs[2:end], barvals))

    return barvals
end

# generate pdf
function generateCPDF(tindx, chosentimeseries, Rs, numSizeClasses)
    # input: tindx = time
    #expanded_data = Vector{Float64}()  # Pre-allocate an empty vector
    barvals = zeros(length(Rs)-1)
    for i in 2:numSizeClasses
        
        nivals = find_z_mean(getfield(chosentimeseries, Symbol("n$i")))  # Extract tracer values
        println(i, nivals[tindx])
        #append!(expanded_data, Rs[i] * nivals[tindx])#ones(round(Int, nivals[tindx])))
        barvals[i-1] = nivals[tindx] * 2*π*Rs[i]^3/(Main.aspect_ratio)
    end
        
    # Perform Kernel Density Estimation of everything apart from size class 1
    #pdfsol = kde(expanded_data)
    #pdfsol = kde((Rs[2:end], barvals))

    return barvals
end

end