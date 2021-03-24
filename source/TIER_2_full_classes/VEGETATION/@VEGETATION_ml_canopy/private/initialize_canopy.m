function [ground] = initialize_canopy(ground)

nleaf = ground.STATVAR.vegetation.mlcanopyinst.nleaf;
ncan = ground.STATVAR.vegetation.mlcanopyinst.ncan;

% ground.STATVAR.vegetation.canopy.lai = ones(1);
% ground.STATVAR.vegetation.canopy.nveg = ones(1);
ground.STATVAR.vegetation.canopy.rd = zeros(1,ncan,nleaf);
% ground.STATVAR.vegetation.canopy.dpai = zeros(1,ncan);

ground.STATVAR.vegetation.canopy.num_exposedvegp = 1;       % Number of non-snow-covered veg points in CLM patch filter
ground.STATVAR.vegetation.canopy.filter_exposedvegp = 1;    % CLM patch filter for non-snow-covered ground.STATVAR.vegetation

% % ground.STATVAR.vegetation.canopy.nsoi = ground.STATVAR.vegetation.soilvar.nsoi;
% ground.STATVAR.vegetation.canopy.nbot = ones(1);
% ground.STATVAR.vegetation.canopy.ntop = ones(1);
% ground.STATVAR.vegetation.canopy.clumpfac = ones(1);
% ground.STATVAR.vegetation.canopy.sumlai = zeros(1,ncan);
% % ground.STATVAR.vegetation.canopy.dlai = zeros(1);
% ground.STATVAR.vegetation.canopy.pai = ones(1);

% --- Define plant canopy
% Set canopy LAI, layer LAI increment, and number of layers
    
%     switch canopy_type
%         case 'dense'
%             lai_inc = 0.1;                                                          % Leaf area index for each layer
%             ground.STATVAR.vegetation.canopy.lai(p) = 6; %used to be 6                             % Leaf area index of canopy (m2/m2)
%         case 'sparse'
%             lai_inc = 0.05;                                                         % Leaf area index for each layer
%             ground.STATVAR.vegetation.canopy.lai(p) = 1;  %used to be 1                            % Leaf area index of canopy (m2/m2)
%     end
%     ground.STATVAR.vegetation.canopy.nveg(p) = round(ground.STATVAR.vegetation.canopy.lai(p) / lai_inc);          % Number of leaf layers in canopy
%     

% % Minimum number of layers for Norman radiation
% 
% if (ground.STATVAR.vegetation.canopy.nveg(p) < 9)
%     ground.STATVAR.vegetation.canopy.nveg(p) = 9;
%     lai_inc = ground.STATVAR.vegetation.canopy.lai(p) / ground.STATVAR.vegetation.canopy.nveg(p);
% end

% % Set array indices for canopy layers
% 
% ground.STATVAR.vegetation.soilvar.nsoi(p) = 1;                               % First layer is soil
% ground.STATVAR.vegetation.canopy.nbot(p) = ground.STATVAR.vegetation.soilvar.nsoi(p) + 1;                 % Bottom leaf layer
% ground.STATVAR.vegetation.canopy.ntop(p) = ground.STATVAR.vegetation.canopy.nbot(p) + ground.STATVAR.vegetation.canopy.nveg(p) - 1;   % Top leaf layer

% % Set LAI of each layer
% 
% for iv = ground.STATVAR.vegetation.canopy.nbot(p):ground.STATVAR.vegetation.canopy.ntop(p)
%     ground.STATVAR.vegetation.canopy.dlai(p,iv) = lai_inc;
% end
% 
% % Cumulative leaf area index (from canopy top) at mid-layer
% 
% for iv = ground.STATVAR.vegetation.canopy.ntop(p): -1: ground.STATVAR.vegetation.canopy.nbot(p)
%     if (iv == ground.STATVAR.vegetation.canopy.ntop(p))
%         ground.STATVAR.vegetation.canopy.sumlai(p,iv) = 0.5 * ground.STATVAR.vegetation.canopy.dlai(p,iv);
%     else
%         ground.STATVAR.vegetation.canopy.sumlai(p,iv) = ground.STATVAR.vegetation.canopy.sumlai(p,iv+1) + ground.STATVAR.vegetation.canopy.dlai(p,iv);
%     end
% end


end