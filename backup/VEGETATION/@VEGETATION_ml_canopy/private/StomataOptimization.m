function [vegetation] = StomataOptimization (vegetation,p,ic,il)
% %DESCRIPTION:
% Photosynthesis and stomatal conductance with optimization

% %ARGUMENTS:
%     implicit none
%     integer, intent(in) :: p              % Patch index for CLM g/l/c/p hierarchy
%     integer, intent(in) :: ic             % Aboveground layer index
%     integer, intent(in) :: il             % Sunlit (1) or shaded (2) leaf index
%     type(mlcanopy_type), intent(inout) :: mlcanopy_inst
%
% %LOCAL VARIABLES:
%     real(r8) :: gs1, gs2                  % Initial guess for gs (mol H2O/m2/s)
%     real(r8) :: check1, check2            % Water-use efficiency and cavitation check for gs1 and gs2

%---------------------------------------------------------------------
% psil        = vegetation.mlcanopyinst.psil   ;  % Leaf water potential (MPa)
% gs          = vegetation.mlcanopyinst.gs     ;  % Leaf stomatal conductance (mol H2O/m2 leaf/s)
% leafwp      = vegetation.mlcanopyinst.leafwp;

% Initialize leaf water potential (for sunlit or shaded leaf) to the layer
% value of the previous time step

vegetation.mlcanopyinst.psil(p,ic,il) = vegetation.mlcanopyinst.lwp(p,ic);

% Low and high initial estimates for gs (mol H2O/m2/s)

gs1 = 0.002;
gs2 = 2;

% Calculate gs

if (vegetation.canopy.dpai(p,ic) > 0.) % leaf layer
    
    % Check for minimum stomatal conductance linked to low light or drought stress
    % based on the water-use efficiency and cavitation checks for gs1 and gs2
    
    [vegetation, check1] = StomataEfficiency (vegetation, p,ic,il, gs1);
    [vegetation, check2] = StomataEfficiency (vegetation, p,ic,il, gs2);
    
    
    if ((check1 * check2) < 0.)
        
        func_name = 'StomataEfficiency';
        tol = 0.004;
        
        % %         % Calculate gs using the function StomataEfficiency to iterate gs
        % %         % to an accuracy of tol (mol H2O/m2/s)
        
        [vegetation, vegetation.mlcanopyinst.gs(p,ic,il)] = brent_root (func_name, vegetation, p, ic, il, gs1, gs2, tol); %gs(p,ic,il) %zbrent?        
        
        % %
        % %         [gs(p,ic,il)] = zbrent (StomataEfficiency, vegetation, gs1, gs2, tol);
        % %
    else
        % %
        % %         % Low light or drought stress. Set gs to minimum conductance
        % %
        vegetation.mlcanopyinst.gs(p,ic,il) = 0.002;
        % %
    end
    
%  vegetation.mlcanopyinst.gs(p,ic,il) = 0.02;
else % non-leaf layer
    
    vegetation.mlcanopyinst.gs(p,ic,il) = 0.;
    
end

% Leaf fluxes and leaf water potential for this gs

[vegetation, vegetation.mlcanopyinst.psil(p,ic,il)] = StomataFluxes (vegetation, p,ic,il, vegetation.mlcanopyinst.gs(p,ic,il), vegetation.mlcanopyinst.psil(p,ic,il));

% vegetation.mlcanopyinst.leafwp      = leafwp ;  % Leaf water potential of canopy layer (MPa)
% vegetation.mlcanopyinst.psil        = psil   ;  % Leaf water potential (MPa)
% vegetation.mlcanopyinst.gs          = gs     ;  % Leaf stomatal conductance (mol H2O/m2 leaf/s)

end

