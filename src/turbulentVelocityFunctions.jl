# calculate the turbulent velocities

module turbulentVelocityFunctions

export get_vₜ_effective_radius_collision_velocity_parameterisation_num_2,
        get_vₜ_collision_velocity_parameterisation_num_2,
        get_vₜ_averaged_radius_collision_velocity_parameterisation_num_2, 
        get_vₜ_effective_radius_collision_velocity_parameterisation_num_13,
        get_vₜ_collision_velocity_parameterisation_num_13,
        get_vₜ_averaged_radius_collision_velocity_parameterisation_num_13,
        get_vₜ_same_radius_averaged_radius_collision_velocity_parameterisation_num_2,
        get_vₜ_same_radius_effective_radius_collision_velocity_parameterisation_num_2,
        get_vₜ_same_radius_collision_velocity_parameterisation_num_2,
        get_vₜ_same_radius_effective_radius_collision_velocity_parameterisation_num_13,
        get_vₜ_same_radius_collision_velocity_parameterisation_num_13,
        get_vₜ_same_radius_averaged_radius_collision_velocity_parameterisation_num_13

# cylindrical approximation different radii
function get_vₜ_effective_radius_collision_velocity_parameterisation_num_2(indx, Rindx, ϵ, ν)
    aspect_ratio = Main.aspect_ratio
    Rs = Main.Rs
    return (3/(2*aspect_ratio))^(1/3) * sqrt(ϵ/(15ν)*(sqrt(π/2) + sqrt(2/π)))*(Rs[indx] + Rs[Rindx]) # cylindrical approximation
end

function get_vₜ_collision_velocity_parameterisation_num_2(indx, Rindx, ϵ, ν)
    Rs = Main.Rs
    return sqrt(ϵ/(15ν)*(sqrt(π/2) + sqrt(2/π)))*(Rs[indx] + Rs[Rindx]) # cylindrical approximation
end

function get_vₜ_averaged_radius_collision_velocity_parameterisation_num_2(indx, Rindx, ϵ, ν)
    effective_areas = Main.effective_areas
    return sqrt(ϵ/(15ν)*(sqrt(π/2) + sqrt(2/π)))*(effective_areas[indx, Rindx]) # cylindrical approximation
end

# spherical approximation different radii
function get_vₜ_effective_radius_collision_velocity_parameterisation_num_13(indx, Rindx, ϵ, ν)
    aspect_ratio = Main.aspect_ratio
    Rs = Main.Rs
    return (3/(2*aspect_ratio))^(1/3) * sqrt(ϵ/(15ν))*(Rs[indx] + Rs[Rindx]) # spherical approximation
end

function get_vₜ_collision_velocity_parameterisation_num_13(indx, Rindx, ϵ, ν)
    Rs = Main.Rs
    return sqrt(ϵ/(15ν))*(Rs[indx] + Rs[Rindx]) # spherical approximation
end

function get_vₜ_averaged_radius_collision_velocity_parameterisation_num_13(indx, Rindx, ϵ, ν)
    effective_areas = Main.effective_areas
    return sqrt(ϵ/(15ν))*(effective_areas[indx, Rindx]) # spherical approximation
end

# cylindrical approximation same radii
function get_vₜ_same_radius_averaged_radius_collision_velocity_parameterisation_num_2(indx, ϵ, ν)
    effective_areas = Main.effective_areas
    return sqrt(ϵ/(15ν)*(sqrt(π/2) + sqrt(2/π)))*(effective_areas[indx, indx]) # cylindrical approximation
end

function get_vₜ_same_radius_effective_radius_collision_velocity_parameterisation_num_2(indx, ϵ, ν)
    aspect_ratio = Main.aspect_ratio
    Rs = Main.Rs
    return (3/(2*aspect_ratio))^(1/3) * sqrt(ϵ/(15ν)*(sqrt(π/2) + sqrt(2/π)))*(2*Rs[indx]) # cylindrical approximation
end

function get_vₜ_same_radius_collision_velocity_parameterisation_num_2(indx, ϵ, ν)
    Rs = Main.Rs
    return sqrt(ϵ/(15ν)*(sqrt(π/2) + sqrt(2/π)))*(2*Rs[indx]) # cylindrical approximation
end

# spherical approximation same radii
function get_vₜ_same_radius_effective_radius_collision_velocity_parameterisation_num_13(indx, ϵ, ν)
    aspect_ratio = Main.aspect_ratio
    Rs = Main.Rs
    return (3/(2*aspect_ratio))^(1/3) * sqrt(ϵ/(15ν))*(2*Rs[indx]) # spherical approximation
end

function get_vₜ_same_radius_collision_velocity_parameterisation_num_13(indx, ϵ, ν)
    Rs = Main.Rs
    return  sqrt(ϵ/(15ν))*(2*Rs[indx]) # spherical approximation
end

function get_vₜ_same_radius_averaged_radius_collision_velocity_parameterisation_num_13(indx, ϵ, ν)
    effective_areas = Main.effective_areas
    return  sqrt(ϵ/(15ν))*(effective_areas[indx, indx]) # spherical approximation
end


end