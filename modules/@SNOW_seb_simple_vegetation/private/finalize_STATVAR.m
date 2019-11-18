function ground = finalize_STATVAR(ground, forcing)

% ground.STATVAR.current_t = 0.0;
ground.STATVAR.execution_t = forcing.PARA.start_time + 0.5;

ground.STATVAR.Lstar = -100;
ground.STATVAR.Qh = 0;
ground.STATVAR.Qe = 0;

% Set up the canopy structure
% [vegetation] = SetUpCanopy();
[ground] = initialize_mlcanopyinst(ground, forcing);
[ground] = initialize_physcon(ground);
[ground] = initialize_params(ground);
[ground] = initialize_pftcon(ground);
[ground] = initialize_leaf(ground);
[ground] = initialize_soil(ground);
[ground] = initialize_atmos(ground);
[ground] = initialize_canopy(ground);
[ground] = initialize_flux(ground);

%build lookup tables here, otherwise the file will be read with each
%iteration step - consumes alot of time
[ground.STATVAR.vegetation.PsiLookup.dtLgridM,...
    ground.STATVAR.vegetation.PsiLookup.zdtgridM,...
    ground.STATVAR.vegetation.PsiLookup.psigridM,...
    ground.STATVAR.vegetation.PsiLookup.dtLgridH,...
    ground.STATVAR.vegetation.PsiLookup.zdtgridH,...
    ground.STATVAR.vegetation.PsiLookup.psigridH] = LookupPsihatINI();
end