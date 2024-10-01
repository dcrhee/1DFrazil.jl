# defines all the constants
module Constants
export Tf, ρₐ, ρₒ, ρᵢ, Cd, cᴾ, Nu, kl, α, Lat

const Tf = -1.6378 # freezing temperature

# densities
const ρₐ = 1.225 # kg/m³ density of air
const ρₒ = 1035 # kg/m³ density of water (might need to change this)
const ρᵢ = 920 # kg/m³ density of ice (might need to change this)

# drag, head capacity
const Cd = 1e-3 # drag coefficient
const cᴾ = 3991.0 # J K⁻¹ kg⁻¹, typical heat capacity for seawater

# Nu
const Nu = 1
const kl = 0.564

const α = 0.31

const Lat = 3.35 * 10^5

end;