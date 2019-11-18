
function [vegetation] = LeafFluxes (vegetation, p, ic, il)
% DESCRIPTION:
% Leaf temperature, energy fluxes, photosynthesis, and stomatal conductance

%   gstyp = vegetation.physcon.gstyp;
%   dpai = vegetation.canopy.dpai;  % Layer plant area index (m2/m2)
%   tleaf = vegetation.mlcanopyinst.tleaf;  % Leaf temperature (K)

if (vegetation.canopy.dpai(p,ic) > 0)   
    % Calculate fluxes for leaf layers using TleafFunc for Ball-Berry style
    % stomatal model or StomataOptimization for water-use efficiency
    % optimization model. These routines calculate leaf temperature and
    % stomatal conductance simultaneously. If leaf_temp_iter = false, the leaf
    % temperature calculation is turned off and the stomatal conductance
    % routines use leaf temperature from the previous sub-timestep.
        
    if (vegetation.params.gstyp <= 1)       
        % Initial estimates for leaf temperature
        t0 = vegetation.mlcanopyinst.tleaf(p,ic,il) - 1 ; % Initial estimate for leaf temperature (K)
        t1 = vegetation.mlcanopyinst.tleaf(p,ic,il) + 1 ; % Initial estimate for leaf temperature (K)
        tol = 0.1;                                        % Accuracy tolerance for tleaf (K)
        
        % Solve for tleaf: Use TleafFunc to iterate leaf temperature, energy fluxes,
        % photosynthesis and stomatal conductance. This temperature is refined to an
        % accuracy of tol. Do not use the returned value (dummy), and instead use
        % the tleaf calculated in the final call to TleafFunc.
        func_name = 'TleafFunc';      % The function name
        vegetation = hybrid_root(func_name, vegetation, p, ic, il, t0, t1, tol);    
        % FORTRAN: dummy = hybrid ('LeafFluxes', p, ic, il, mlcanopy_inst, TleafFunc, t0, t1, tol)  
        
    elseif (vegetation.params.gstyp == 2)       
        % Iterate leaf temperature, flux calculations, and stomatal conductance
        % using water-use efficiency optimization and cavitation check       
        vegetation = StomataOptimization(vegetation,p,ic,il);        
    end 
    
else % non-leaf layer   
     % Zero out fluxes    
    if (vegetation.params.gstyp <= 1)
        vegetation = TleafFunc(vegetation,p,ic,il,vegetation.mlcanopyinst.tair(p,ic));
    elseif (vegetation.params.gstyp == 2)
        vegetation = StomataOptimization(vegetation,p,ic,il);
    end    
end

%   vegetation.mlcanopyinst.t1 = t1;
end
  
