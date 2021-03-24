function [vegetation] = PlantResistance (vegetation) %, p, ic
%num_exposedvegp, filter_exposedvegp, mlcanopy_inst, shr_kind, ic, p, f)
%
% %DESCRIPTION:
% Calculate whole-plant leaf-specific conductance (soil-to-leaf)

% %LOCAL VARIABLES:
%     integer  :: f                          % Filter index
% p  = vegetation.mlcanopyinst.p;                       % Patch index for CLM g/l/c/p hierarchy
% ic = 1;                       % Aboveground layer index
%     real(r8) :: rplant                     % Aboveground plant hydraulic resistance (MPa.s.m2/mmol H2O)

% ncan = vegetation.mlcanopyinst.ncan;
%---------------------------------------------------------------------
%Input
% gplant % Stem (xylem-to-leaf) hydraulic conductance (mmol H2O/m2 leaf area/s/MPa)
% dpai % Layer plant area index (m2/m2)
% zs % Canopy height for scalar concentration and source (m)
% rsoil % Soil hydraulic resistance (MPa.s.m2/mmol H2O)

for f = 1:vegetation.canopy.num_exposedvegp
    p = vegetation.canopy.filter_exposedvegp(f);
    for ic = vegetation.canopy.nbot(p):vegetation.mlcanopyinst.ncan(p) %Fortran 1:ncan
        
        if (vegetation.canopy.dpai(p,ic) > 0) % leaf layer
            
            % Aboveground plant stem resistance, xylem-to-leaf (MPa.s.m2/mmol H2O)
            
            %rplant = zs(p,ic) / gplant(p)       % gplant is conductivity (mmol/m/s/MPa)
            rplant = 1.  ./ vegetation.leaf.gplant(p);         % gplant is conductance (mmol/m2/s/MPa)
            
            % Leaf specific conductance, soil-to-leaf (mmol H2O/m2/s/MPa)
            
            %rsoil calculated in SoilResistance
            vegetation.mlcanopyinst.lsc(p,ic) = 1.  ./ (vegetation.mlcanopyinst.rsoil(p) + rplant);
            
        else % non-leaf layer
            
            vegetation.mlcanopyinst.lsc(p,ic) = 0.;
            
        end
        
    end
end
end
