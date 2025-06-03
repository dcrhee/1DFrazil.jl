module encounterFrequencyFunctions

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

function get_Fenc_effective_radius_collision_velocity_parameterisation_num_3(indx, Rindx, nRindx)
    aspect_ratio = Main.aspect_ratio
    Rs = Main.Rs
    vcoll_matrix = Main.vcoll_matrix
    if indx == Rindx
        Fenc = (3/(2*aspect_ratio))^(2/3) * π*(Rs[indx] + Rs[Rindx])^2 * nRindx * vcoll_matrix[indx, Rindx]
    else
        Fenc = (3/(2*aspect_ratio))^(2/3) * 2*π*(Rs[indx] + Rs[Rindx])^2 * nRindx * vcoll_matrix[indx, Rindx]
    end
    return Fenc
end

function get_Fenc_average_radius_collision_velocity_parameterisation_num_3(indx, Rindx, nRindx)
    vcoll_matrix = Main.vcoll_matrix
    effective_areas = Main.effective_areas
    if indx == Rindx
        Fenc =  π*(effective_areas[indx, Rindx])^2 * nRindx * vcoll_matrix[indx, Rindx]
    else
        Fenc = 2*π*(effective_areas[indx, Rindx])^2 * nRindx * vcoll_matrix[indx, Rindx]
    end
    return Fenc
end

function get_Fenc_collision_velocity_parameterisation_num_3(indx, Rindx, nRindx)
    Rs = Main.Rs
    vcoll_matrix = Main.vcoll_matrix
    if indx == Rindx
        Fenc =  π*(Rs[indx] + Rs[Rindx])^2 * nRindx * vcoll_matrix[indx, Rindx]
    else
        Fenc =  2*π*(Rs[indx] + Rs[Rindx])^2 * nRindx * vcoll_matrix[indx, Rindx]
    end
    return Fenc
end


function get_Fenc_average_radius_collision_velocity_parameterisation_num_12(indx, Rindx, nRindx)
    vcoll_matrix = Main.vcoll_matrix
    effective_areas = Main.effective_areas
    if indx == Rindx
        Fenc =  π * (effective_areas[indx, Rindx])^2 * nRindx/2 * vcoll_matrix[indx, Rindx]
    else
        Fenc = π * (effective_areas[Rindx, Rindx])^2 * nRindx * vcoll_matrix[indx, Rindx]
    end
    return Fenc
end

function get_Fenc_effective_radius_collision_velocity_parameterisation_num_12(indx, Rindx, nRindx)
    aspect_ratio = Main.aspect_ratio
    Rs = Main.Rs
    vcoll_matrix = Main.vcoll_matrix
    if indx == Rindx
        Fenc = (3/(2*aspect_ratio))^(2/3) * π*(Rs[indx] + Rs[Rindx])^2 * nRindx/2 * vcoll_matrix[indx, Rindx]
    else
        Fenc = (3/(2*aspect_ratio))^(2/3) * π*(Rs[indx] + Rs[Rindx])^2 * nRindx * vcoll_matrix[indx, Rindx]
    end
    return Fenc
end

function get_Fenc_collision_velocity_parameterisation_num_12(indx, Rindx, nRindx)
    Rs = Main.Rs
    vcoll_matrix = Main.vcoll_matrix
    if indx == Rindx
        Fenc = π*(Rs[indx] + Rs[Rindx])^2 * nRindx/2 * vcoll_matrix[indx, Rindx]
    else
        Fenc = π*(Rs[indx] + Rs[Rindx])^2 * nRindx * vcoll_matrix[indx, Rindx]
    end
    return Fenc
end

function get_Fenc_p1_average_radius_collision_velocity_parameterisation_num_3(indx, ntot)
    vcoll_matrix = Main.vcoll_matrix
    effective_areas = Main.effective_areas
    #Fenc =  (3/(2*aspect_ratio))^(2/3) * 2*π*(Rs[indx])^2 * vcoll * ntot
    return  2*π*(effective_areas[indx, indx]/2)^2 * vcoll_matrix[indx] * ntot
end

function get_Fenc_p1_effective_radius_collision_velocity_parameterisation_num_3(indx, ntot)
    aspect_ratio = Main.aspect_ratio
    Rs = Main.Rs
    vcoll_matrix = Main.vcoll_matrix
    #Fenc =  (3/(2*aspect_ratio))^(2/3) * 2*π*(Rs[indx])^2 * vcoll * ntot
    return (3/(2*aspect_ratio))^(2/3) * 2*π*(Rs[indx])^2 * vcoll_matrix[indx] * ntot
end

function get_Fenc_p1_collision_velocity_parameterisation_num_3(indx, ntot)
    Rs = Main.Rs
    vcoll_matrix = Main.vcoll_matrix
    #Fenc = 2*π*(Rs[indx])^2 * vcoll * ntot
    return 2*π*(Rs[indx])^2 * vcoll_matrix[indx] * ntot
end

function get_Fenc_p1_average_radius_collision_velocity_parameterisation_num_12(indx, ntot)
    vcoll_matrix = Main.vcoll_matrix
    effective_areas = Main.effective_areas
    #Fenc = (3/(2*aspect_ratio))^(2/3) *π*(Rs[indx])^2 * vcoll * ntot
    return π * (effective_areas[indx, indx]/2)^2  * vcoll_matrix[indx] * ntot
end

function get_Fenc_p1_effective_radius_collision_velocity_parameterisation_num_12(indx, ntot)
    aspect_ratio = Main.aspect_ratio
    Rs = Main.Rs
    vcoll_matrix = Main.vcoll_matrix
    #Fenc = (3/(2*aspect_ratio))^(2/3) *π*(Rs[indx])^2 * vcoll * ntot
    return (3/(2*aspect_ratio))^(2/3) *π*(Rs[indx])^2 * vcoll_matrix[indx] * ntot
end

function get_Fenc_p1_collision_velocity_parameterisation_num_12(indx, ntot)
    Rs = Main.Rs
    vcoll_matrix = Main.vcoll_matrix
    #Fenc = π*(Rs[indx])^2 * vcoll * ntot
    return π*(Rs[indx])^2 * vcoll_matrix[indx] * ntot
end

end