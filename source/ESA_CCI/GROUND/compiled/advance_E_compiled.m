
function energy = advance_E_compiled(energy, d_energy, timestep)
            
energy = energy + d_energy .* timestep;