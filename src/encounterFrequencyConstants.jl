module encounterFrequencyConstants

export get_Fenc_effective_radius_collision_velocity_parameterisation_num_3,
        get_Fenc_average_radius_collision_velocity_parameterisation_num_3,
        get_Fenc_collision_velocity_parameterisation_num_3,

        get_Fenc_average_radius_collision_velocity_parameterisation_num_12,
        get_Fenc_effective_radius_collision_velocity_parameterisation_num_12,
        get_Fenc_collision_velocity_parameterisation_num_12,

        get_Fenc_p1_average_radius_collision_velocity_parameterisation_num_3,
        get_Fenc_p1_effective_radius_collision_velocity_parameterisation_num_3,
        get_Fenc_p1_collision_velocity_parameterisation_num_3,

        get_Fenc_p1_average_radius_collision_velocity_parameterisation_num_12,
        get_Fenc_p1_effective_radius_collision_velocity_parameterisation_num_12,
        get_Fenc_p1_collision_velocity_parameterisation_num_12

function get_Fenc_effective_radius_collision_velocity_parameterisation_num_3(indx, Rindx)
    aspect_ratio = Main.aspect_ratio
    Rs = Main.Rs
    vcoll_matrix = Main.vcoll_matrix
    if indx == Rindx
        Fenc = (3/(2*aspect_ratio))^(2/3) * π*(Rs[indx] + Rs[Rindx])^2 * vcoll_matrix[indx, Rindx]
    else
        Fenc = (3/(2*aspect_ratio))^(2/3) * 2*π*(Rs[indx] + Rs[Rindx])^2 * vcoll_matrix[indx, Rindx]
    end
    return Fenc
end

function get_Fenc_average_radius_collision_velocity_parameterisation_num_3(indx, Rindx)
    Rs = Main.Rs
    vcoll_matrix = Main.vcoll_matrix
    radii_efficiency = Main.radii_efficiency
    if indx == Rindx
        Fenc =  π*(Rs[indx] + Rs[Rindx])^2*(radii_efficiency[indx, Rindx]) * vcoll_matrix[indx, Rindx]
    else
        Fenc = 2*π*(Rs[indx] + Rs[Rindx])^2*(radii_efficiency[indx, Rindx]) * vcoll_matrix[indx, Rindx]
    end
    return Fenc
end

function get_Fenc_collision_velocity_parameterisation_num_3(indx, Rindx)
    Rs = Main.Rs
    vcoll_matrix = Main.vcoll_matrix
    if indx == Rindx
        Fenc =  π*(Rs[indx] + Rs[Rindx])^2 * vcoll_matrix[indx, Rindx]
    else
        Fenc =  2*π*(Rs[indx] + Rs[Rindx])^2 * vcoll_matrix[indx, Rindx]
    end
    return Fenc
end


function get_Fenc_average_radius_collision_velocity_parameterisation_num_12(indx, Rindx)
    Rs = Main.Rs
    vcoll_matrix = Main.vcoll_matrix
    radii_efficiency = Main.radii_efficiency
    if indx == Rindx
        Fenc =  π * (Rs[indx] + Rs[Rindx])^2/2 * (radii_efficiency[indx, Rindx]) * vcoll_matrix[indx, Rindx]
    else
        Fenc = π * (Rs[indx] + Rs[Rindx])^2 * (radii_efficiency[Rindx, Rindx]) * vcoll_matrix[indx, Rindx]
    end
    return Fenc
end

function get_Fenc_effective_radius_collision_velocity_parameterisation_num_12(indx, Rindx)
    aspect_ratio = Main.aspect_ratio
    Rs = Main.Rs
    vcoll_matrix = Main.vcoll_matrix
    if indx == Rindx
        Fenc = (3/(2*aspect_ratio))^(2/3) * π*(Rs[indx] + Rs[Rindx])^2/2 * vcoll_matrix[indx, Rindx]
    else
        Fenc = (3/(2*aspect_ratio))^(2/3) * π*(Rs[indx] + Rs[Rindx])^2 * vcoll_matrix[indx, Rindx]
    end
    return Fenc
end

function get_Fenc_collision_velocity_parameterisation_num_12(indx, Rindx)
    Rs = Main.Rs
    vcoll_matrix = Main.vcoll_matrix
    if indx == Rindx
        Fenc = π*(Rs[indx] + Rs[Rindx])^2/2 * vcoll_matrix[indx, Rindx]
    else
        Fenc = π*(Rs[indx] + Rs[Rindx])^2 * vcoll_matrix[indx, Rindx]
    end
    return Fenc
end

function get_Fenc_p1_average_radius_collision_velocity_parameterisation_num_3(indx)
    Rs = Main.Rs
    vcoll_matrix = Main.vcoll_matrix
    radii_efficiency = Main.radii_efficiency
    #Fenc =  (3/(2*aspect_ratio))^(2/3) * 2*π*(Rs[indx])^2 * vcoll * ntot
    return  2*π*(Rs[indx])^2 * (radii_efficiency[indx, indx]) * vcoll_matrix[indx]
end

function get_Fenc_p1_effective_radius_collision_velocity_parameterisation_num_3(indx)
    aspect_ratio = Main.aspect_ratio
    Rs = Main.Rs
    vcoll_matrix = Main.vcoll_matrix
    #Fenc =  (3/(2*aspect_ratio))^(2/3) * 2*π*(Rs[indx])^2 * vcoll * ntot
    return (3/(2*aspect_ratio))^(2/3) * 2*π*(Rs[indx])^2 * vcoll_matrix[indx]
end

function get_Fenc_p1_collision_velocity_parameterisation_num_3(indx)
    Rs = Main.Rs
    vcoll_matrix = Main.vcoll_matrix
    #Fenc = 2*π*(Rs[indx])^2 * vcoll * ntot
    return 2*π*(Rs[indx])^2 * vcoll_matrix[indx]
end

function get_Fenc_p1_average_radius_collision_velocity_parameterisation_num_12(indx)
    Rs = Main.Rs
    vcoll_matrix = Main.vcoll_matrix
    radii_efficiency = Main.radii_efficiency
    #Fenc = (3/(2*aspect_ratio))^(2/3) *π*(Rs[indx])^2 * vcoll * ntot
    return π*(Rs[indx])^2 * vcoll_matrix[indx] *  (radii_efficiency[indx, indx])
end

function get_Fenc_p1_effective_radius_collision_velocity_parameterisation_num_12(indx)
    aspect_ratio = Main.aspect_ratio
    Rs = Main.Rs
    vcoll_matrix = Main.vcoll_matrix
    #Fenc = (3/(2*aspect_ratio))^(2/3) *π*(Rs[indx])^2 * vcoll * ntot
    return (3/(2*aspect_ratio))^(2/3) *π*(Rs[indx])^2 * vcoll_matrix[indx]
end

function get_Fenc_p1_collision_velocity_parameterisation_num_12(indx)
    Rs = Main.Rs
    vcoll_matrix = Main.vcoll_matrix
    #Fenc = π*(Rs[indx])^2 * vcoll * ntot
    return π*(Rs[indx])^2 * vcoll_matrix[indx]
end

end