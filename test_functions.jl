function find_density(T, S, ρ₀ = 1027.0, β = 0.0078, α = 1.67*10^(-4))
    ρ = ρ₀ * (1.0256550500000001 - α*T + β*S)
    return ρ
end

function temperature_constant(Tₒ, Sₒ, Rᵢ, H = 0.0004) # inputs are T and S at time = 0 as we assume the density stays constant
    # constant in front of each concentration
    ρ = find_density(Tₒ, Sₒ)
    Tconstᵢ = kl*Nu/(ρ*cᴾ) *2π * Rᵢ *  1/(0.9002 - 0.2634*log(H/(2*Rᵢ)))
    return Tconstᵢ
end

function S_const(T₀) # inputs are T and S at time = 0 as we assume the density stays constant
    # constant in front of each concentration
    Sconst = -(1-αₛ)*(Tf - T₀)cᴾ/Lat
    return Sconst
end

function T_const(Tₒ, Sₒ, n₁ₒ, n₂ₒ, n₃ₒ, R₁, R₂, R₃) # input the concentrations at t = 0
    Tconst₁ = temperature_constant(Tₒ, Sₒ, R₁)
    Tconst₂ = temperature_constant(Tₒ, Sₒ, R₂)
    Tconst₃ = temperature_constant(Tₒ, Sₒ, R₃)
    return Tconst₁ * n₁ₒ + Tconst₂ * n₂ₒ + Tconst₃ * n₃ₒ
end

function G_const(Rᵢ, H = 0.0004)
    Gᵢ = kl*Nu/(ρᵢ*Lat) *2π * Rᵢ *  1/(0.9002 - 0.2634*log(H/(2*Rᵢ)))
    return Gᵢ
end

function find_Vi(R, H = 0.0004)
    return π*R^2*H
end