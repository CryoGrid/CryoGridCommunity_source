% function E = advance_E_compiled(E, timestep, dE_dt)
% E = E + timestep .* dE_dt;

function energy = advance_E_compiled(energy, d_energy, timestep)
            
energy = energy + d_energy .* timestep;