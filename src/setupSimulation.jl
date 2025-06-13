module setupSimulation
# progress simulation functions
using Printf: @sprintf
using Oceananigans
export progress, enforce_nonnegative_tracer

function progress(simulation)
    u, v, w = simulation.model.velocities
    T = simulation.model.tracers.T
    n1 = simulation.model.tracers.n1
    n2 = simulation.model.tracers.n2
    n3 = simulation.model.tracers.n3
    n4 = simulation.model.tracers.n4
    n5 = simulation.model.tracers.n5
    n6 = simulation.model.tracers.n6
    n7 = simulation.model.tracers.n7
    n8 = simulation.model.tracers.n8
    n9 = simulation.model.tracers.n9
    n10 = simulation.model.tracers.n10
    n50 = simulation.model.tracers.n50 
    n100 = simulation.model.tracers.n100
    n200 = simulation.model.tracers.n200

    # Print a progress message
    #msg = @sprintf("i: %04d, t: %s, Δt: %s, umax = (%.1e, %.1e, %.1e) ms⁻¹, wall time: %s\n",
    msg = @sprintf("i: %04d, t: %s, Δt: %s, umax = (%.1e, %.1e, %.1e) ms⁻¹, n1 = %.1e, n2 = %.1e, n3 = %.1e, n4 = %.1e, n5 = %.1e, n6 = %.1e, n7 = %.1e, n8 = %.1e, n9 = %.1e, n10 = %.1e, n50 = %.1e, n100 = %.1e, n200 = %.1e, Tmin = %.5f, Tmax = %.5f, wall time: %s\n",
    iteration(simulation),
    prettytime(time(simulation)),
    prettytime(simulation.Δt),
    maximum(abs, u), maximum(abs, v), maximum(abs, w),
    minimum(n1), minimum(n2), minimum(n3), minimum(n4), minimum(n5), minimum(n6), minimum(n7), minimum(n8), minimum(n9), minimum(n10), minimum(n50), minimum(n100), minimum(n200), 
    #maximum(n1), maximum(n2), maximum(n3), maximum(n4), maximum(n5), maximum(n6), maximum(n7), maximum(n8), maximum(n9), maximum(n10), 
    minimum(T), maximum(T),
    prettytime(simulation.run_wall_time))

    @info msg

    return nothing
end

function enforce_nonnegative_tracer(simulation)
    for n in 1:Main.numSizeClasses
        tracer_data = getproperty(simulation.model.tracers, Symbol("n$n")).data
        #@inbounds tracer_data .= max.(tracer_data, 1e-50)
        @inbounds tracer_data .= max.(tracer_data, 0)
    end
    return nothing
end

end