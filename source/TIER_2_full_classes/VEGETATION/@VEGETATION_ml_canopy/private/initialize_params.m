
function [ground] = initialize_params(ground)

%clm_arguments structure/container
ground.STATVAR.vegetation.params.numrad = 2;         % Number of wavebands (visible, near-infrared)
ground.STATVAR.vegetation.params.vis = 1;            % Array index for visible waveband
ground.STATVAR.vegetation.params.nir = 2;            % Array index for near-infrared waveband
ground.STATVAR.vegetation.params.sun = 1;            % Array index for sunlit leaf
ground.STATVAR.vegetation.params.sha = 2;            % Array index for shaded leaf
ground.STATVAR.vegetation.params.npts = 1;           % Number of grid points to process

ground.STATVAR.vegetation.params.gstyp = 2; %1;                    % Stomatal conductance: 0 = Medlyn model. 1 = Ball-Berry model. 2 = WUE optimization

end
