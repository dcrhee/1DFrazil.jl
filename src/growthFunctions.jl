#growthFunctions

module growthFunctions

export find_growth_rate, find_density

function find_growth_rate(T, Rᵢ, H)
    #G = kl*Nu/(ρᵢ*Lat) * (Tf - T) *2π * Rᵢ *  1/(0.9002 - 0.2634*log(H/(2*Rᵢ)))
    G = kl*Nu/(ρᵢ*Lat) * (Tf - T) * 2π * H
    return G
end

function find_density(T, S)
    ρ = ρₒ * (1 + 7.86 * 1e-4 * (S - 34.5) - 3.87 * 1e-5 * (T + 2) )
    return ρ
end

end
