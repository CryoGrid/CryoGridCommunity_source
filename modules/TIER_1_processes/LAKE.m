%========================================================================
% CryoGrid TIER1 library class for functions related to water bodies 
% S. Westermann, October 2020
%========================================================================

classdef LAKE < BASE

    methods
        
        %move melted grid cells below ice cover (frozen grid cells)
        function ground = move_ice_up(ground)
            fully_melted = (ground.STATVAR.energy > 0);
            ground.STATVAR.energy = reorganize_cells_frozen_melted(ground, ground.STATVAR.energy, fully_melted);
            ground.STATVAR.waterIce = reorganize_cells_frozen_melted(ground, ground.STATVAR.waterIce, fully_melted);
            ground.STATVAR.mineral = reorganize_cells_frozen_melted(ground, ground.STATVAR.mineral, fully_melted);
            ground.STATVAR.organic = reorganize_cells_frozen_melted(ground, ground.STATVAR.organic, fully_melted);
            %ground.STATVAR.air = reorganize_cells_frozen_melted(ground, ground.STATVAR.air, fully_melted);
            ground.STATVAR.layerThick = reorganize_cells_frozen_melted(ground, ground.STATVAR.layerThick, fully_melted);
        end
        
        function reorganized = reorganize_cells_frozen_melted(ground, variable, fully_melted)
            melted_cells = variable(fully_melted); 
            frozen_cells = variable(~fully_melted);
            reorganized = [frozen_cells; melted_cells];
        end
        
        %stratifies the water column, move 1/2 cell up per timestep when water density is lower
        function ground = stratify(ground)  
            if size(ground.STATVAR.energy,1) > 1
                density_water = water_density(ground);
                swap = double(ground.STATVAR.energy(1:end-1,1) >=0 & ground.STATVAR.energy(2:end,1) >= 0 & density_water(2:end,1) < density_water(1:end-1,1));
                
                %only energy is moved, water not moved physically, change if there are solutes, etc.
                energy_down = swap .* 0.5 .* min(ground.STATVAR.waterIce(1:end-1,1), ground.STATVAR.waterIce(2:end,1)) ./ ground.STATVAR.waterIce(1:end-1,1) .* ground.STATVAR.energy(1:end-1,1);
                energy_up = swap .* 0.5 .* min(ground.STATVAR.waterIce(1:end-1,1), ground.STATVAR.waterIce(2:end,1)) ./ ground.STATVAR.waterIce(2:end,1) .* ground.STATVAR.energy(2:end,1);
                
                ground.STATVAR.energy(1:end-1,1) = ground.STATVAR.energy(1:end-1,1) - energy_down + energy_up;
                ground.STATVAR.energy(2:end,1) = ground.STATVAR.energy(2:end,1)  + energy_down - energy_up;
            end
        end
        
        %density of water
        function density_water = water_density(ground)
            T = max(0, ground.STATVAR.T);
            %density_water = (999.83952 + 16.945176 .* T - 7.9870401e-3 .* T.^2 - 46.170461e-6 .* T.^3 + 105.56302e-9 .* T.^4 - 280.54253e-12 .* T.^5) ./ (1 + 16.897850e-3 .* T);
            %from Kell 1975, using simplified version, works well 0-30 degreeC 
            density_water = (999.83952 + 16.945176 .* T - 7.9870401e-3 .* T.^2 - 46.170461e-6 .* T.^3) ./ (1 + 16.897850e-3 .* T);
        end
        
    end
end

