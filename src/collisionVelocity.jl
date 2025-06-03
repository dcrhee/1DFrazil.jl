# collision functions
module collisionVelocity

using SpecialFunctions

export find_vcoll

function find_vcoll(vᵣ, vₜ, collision_velocity_parameterisation_num)
    if collision_velocity_parameterisation_num == 1
        vcoll = sqrt(vᵣ^2 + vₜ^2)
    elseif collision_velocity_parameterisation_num == 2
        vcoll = sqrt(vᵣ^2 + (π+2)^2/(2π)*vₜ^2)
    else
        if vᵣ == 0
            vcoll = sqrt(2/π)*vₜ
        else
            vcoll = sqrt(2/π)*vₜ*(1/2*exp(-vᵣ^2/(2*vₜ^2)) + sqrt(π)/(2*sqrt(2))*(vᵣ/vₜ + vₜ/vᵣ)*erf(vᵣ/(sqrt(2)*vₜ)))
        end
    end
    return vcoll
end

end