module riseVelocityFunctions

export find_steady_velocity

function find_steady_velocity(indx)
    vgs = Main.vgs
    #uᵢ = 30*Rᵢ^(1.2)
    uᵢ = vgs[indx]
    return uᵢ
end

end