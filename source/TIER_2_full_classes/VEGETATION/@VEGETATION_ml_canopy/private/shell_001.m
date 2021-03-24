clear all 
close all

% Set up the canopy structure
vegetation=[];
[vegetation] = SetUpCanopy(vegetation);

[FORCING] = initialize_forcing(21309);

% Set up the canopy structure
% [vegetation] = SetUpCanopy();
[vegetation] = initialize_mlcanopyinst(vegetation, FORCING);
[vegetation] = initialize_physcon(vegetation);
[vegetation] = initialize_params(vegetation);
[vegetation] = initialize_pftcon(vegetation);
[vegetation] = initialize_leaf(vegetation);
[vegetation] = initialize_soil(vegetation);
[vegetation] = initialize_atmos(vegetation);
[vegetation] = initialize_canopy(vegetation);
[vegetation] = initialize_flux(vegetation);

%build lookup tables here, otherwise the file will be read with each
%iteration step - consumes alot of time
[vegetation.PsiLookup.dtLgridM,...
    vegetation.PsiLookup.zdtgridM,...
    vegetation.PsiLookup.psigridM,...
    vegetation.PsiLookup.dtLgridH,...
    vegetation.PsiLookup.zdtgridH,...
    vegetation.PsiLookup.psigridH] = LookupPsihatINI();

tic()
% counter=1;


%summer
%FORCING.data.t_span(21160) - FORCING.data.t_span(23375) %'01-Jun-1981' - '31-Aug-1981 23:00:00'

for i=21309:21350 %50 %02 % 22925 %22993:23030 %22882:22920 %48897:48936 %22870:22910 %17.534-27.07ï¿½C %22952:22967 %
    
    [FORCING] = initialize_forcing(i);
    
    % Set up required forcing structure
    [vegetation] = SetUpForcing(vegetation, FORCING);
    
    % Run the Bonan Model
    profile on
    %     tic
    [vegetation, p, ic, il] = CanopyFluxesMultilayer(vegetation); %,counter);
%     counter = counter + 1;
    %     toc
    profile off
    
    height = vegetation.mlcanopyinst.zw./vegetation.mlcanopyinst.ztop;
    
    figure (1)
    plot([vegetation.mlcanopyinst.tleaf(:,:,1)'-273.15 vegetation.mlcanopyinst.tair(:,:)'-273.15 vegetation.mlcanopyinst.tveg(:,:,1)'-273.15])
    axis([1 25 10 30])
    % set(gcf,'position',[300 50 400 250])
    % set(gcf,'position',[360 200 600 500])
    title('Temperature of a diurnal cycle')
    xlabel('Canopy layer')
    ylabel('Temperature (°C)')
    pause(0.1)
    hold on
    plot(vegetation.mlcanopyinst.tg'-273.15,'*')
    
    figure (2) %blue red yellow purple
    plot(vegetation.mlcanopyinst.st_prof')
    axis([1 25 -40 100])
    % set(gcf,'position',[300 200 600 450])
    % set(gcf,'position',[360 800 600 500])
    title('Storage heat flux')
    hold on
    pause(0.1)
    
    figure (3)
    plot(vegetation.mlcanopyinst.sh_prof')
    axis([1 25 -40 100])
    % set(gcf,'position',[980 200 600 500])
    title('Sensibel heat fluxes of a diurnal cycle')
    xlabel('Canopy layer')
    ylabel('Flux (W/m^2)')
    hold on
    pause(0.1)
    
    figure (4)
    plot(vegetation.mlcanopyinst.lh_prof')
    axis([1 25 -40 100])
    % set(gcf,'position',[980 800 600 500])
    title('Latent heat fluxes of a diurnal cycle')
    xlabel('Canopy layer')
    ylabel('Flux (W/m^2)')
    hold on
    pause(0.1)
    
    figure (5)
    plot (vegetation.mlcanopyinst.rn_prof')
    axis([1 25 -40 150])
    % set(gcf,'position',[1600 200 600 500])
    title('Net radiation of a diurnal cycle')
    xlabel('Canopy layer')
    ylabel('Flux (W/m^2)')
    hold on
    pause(0.1)
    plot(vegetation.flux.swsoi','*')
    plot(vegetation.flux.irsoi','x')
    
    figure (10)
    plot(vegetation.mlcanopyinst.ustar', '*')
    title('Friction velocity u*')
    hold on
    % set(gcf,'position',[2200 800 600 500])
    pause(0.1)
    
    
end

% bowen2=vegetation.mlcanopyinst.shflx/vegetation.mlcanopyinst.lhflx;
% bow_leaf = vegetation.mlcanopyinst.shleaf./vegetation.mlcanopyinst.lhleaf;
% vegetation.mlcanopyinst.sh_prof./vegetation.mlcanopyinst.lh_prof;
% disp(bowen2);

figure (6)
% wind profile
plot(vegetation.mlcanopyinst.wind,height);
% set(gcf,'position',[1600 800 600 500])
title('Wind profile (last timestep)')
xlabel('Windspeed (m/s^2)')
ylabel('Fraction of canopy height')
hold on
% plot(vegetation.canopy.dlai,height)


figure (7)
% vapor pressure deficit profile
plot(vegetation.mlcanopyinst.vpd(:,:,2),height);
% set(gcf,'position',[2220 200 600 500])
title('Vapor pressure deficit profile (last timestep)')


figure (8)
% canopy temperature profile (last timestep)
plot(vegetation.mlcanopyinst.tair(:,:,1)-273.15,height)
hold
% plot(vegetation.mlcanopyinst.tleaf(:,:,1)-273.15,height)
% plot(vegetation.mlcanopyinst.tleaf(:,:,2)-273.15,height)
plot(vegetation.mlcanopyinst.tveg(:,:,2)-273.15,height)
% set(gcf,'position',[2220 800 600 500])
title('Temperature profile (last timestep)')
xlabel('Temperature (°C)')
ylabel('Fraction of canopy height')


% canopy photosynthesis profile (last timestep)
figure (9)
plot(vegetation.mlcanopyinst.an(:,:,1),height)
% set(gcf,'position',[2840 200 600 500])
title('Leaf photosynthesis rate profile (last timestep)')
hold
plot(vegetation.mlcanopyinst.an(:,:,2),height)
xlabel('Photosynthesis rate (umol CO2/m^2 leaf/s)')
ylabel('Fraction of canopy height')

figure (11)
% canopy shaded / sunny fraction
plot(vegetation.flux.fracsun,height)
title('Fraction of sunny and shaded leaves')
hold
plot(vegetation.flux.fracsha,height)
xlabel('Fraction sunny (blue) vs. shaded (red)')
ylabel('Fraction of canopy height')

% Aerodynamic conductance (last timestep)
% plot(vegetation.mlcanopyinst.ga_prof,height)

toc()

% % %graph function dependencies
% % %open /doc/index.html with browser
% % cd ..
% % addpath m2html
% % system('rm doc')
% % m2html('mfiles','20190522_Vegetation_clm_5_current', 'htmldir','doc', 'recursive','on', 'global','on','template','frame', 'index','menu', 'graph','on');

