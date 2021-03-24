function vegetation = set_up_forcing(vegetation, forcing)

ivis = vegetation.params.vis; % Array index for visible waveband
inir = vegetation.params.nir; % Array index for near-infrared waveband
isun = vegetation.params.sun; % Array index for sunlit leaf
isha = vegetation.params.sha; % Array index for shaded leaf

% vegetation.atmos.solar_zenith = 30*pi/180; %angle from earth axes
% vegetation.sun.azimuth = 90 - El;

for f = 1:vegetation.canopy.num_exposedvegp
    p = vegetation.canopy.filter_exposedvegp(f);
    
    
    [Az, El] = SolarAzEl(forcing.TEMP.t,63.18946,118.19596,121); %52.5243700,13.4105300,43
    El = max(El, 5); %ML: min sun elevation fix
    vegetation.sun.solar_az(p) = Az*pi/180;  % Az Azimuth location of the sun (deg)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    vegetation.sun.solar_zenit(p) = (90-El) * pi / 180;
    
    % Atmospheric forcing: CLM grid cell (g) variables -> patch (p) variables
    
    vegetation.mlcanopyinst.zref(p) = vegetation.mlcanopyinst.zref;  %forc_hgt(g);
    vegetation.mlcanopyinst.uref(p) = forcing.TEMP.wind; %ground.STATVAR.vegetation.mlcanopyinst.uref;  %max (1.0 , sqrt(forc_u(g)*forc_u(g)+forc_v(g)*forc_v(g)));
    
    %SEBAS: shouldn't this be distributed, like below?
    vegetation.atmos.swskyb(p,ivis) = forcing.TEMP.Sin_dir .* 0.7;            % Direct beam solar radiation for visible waveband (W/m2)
    vegetation.atmos.swskyd(p,ivis) = forcing.TEMP.Sin_dif .* 0.7;            % Diffuse solar radiation for visible waveband (W/m2)
    vegetation.atmos.swskyb(p,inir) = forcing.TEMP.Sin_dir .* 0.3;            % Direct beam solar radiation for near-infrared waveband (W/m2)
    vegetation.atmos.swskyd(p,inir) = forcing.TEMP.Sin_dif .* 0.3;            % Diffuse solar radiation for near-infrared waveband (W/m2)   
    
%        vegetation.atmos.swskyb(p,ivis) = forcing.TEMP.Sin*0.9;            % Direct beam solar radiation for visible waveband (W/m2)
%        vegetation.atmos.swskyd(p,ivis) = forcing.TEMP.Sin*0.1;            % Diffuse solar radiation for visible waveband (W/m2)
%        vegetation.atmos.swskyb(p,inir) = forcing.TEMP.Sin*0.9;            % Direct beam solar radiation for near-infrared waveband (W/m2)
%        vegetation.atmos.swskyd(p,inir) = forcing.TEMP.Sin*0.1;            % Diffuse solar radiation for near-infrared waveband (W/m2)
    
    % Atmospheric forcing: CLM column (c) variables -> patch (p) variables
    
    vegetation.mlcanopyinst.tref(p) = forcing.TEMP.Tair+273.15; %ground.STATVAR.vegetation.mlcanopyinst.tref; %forc_t(c);
    vegetation.mlcanopyinst.qref(p) = forcing.TEMP.q; %ground.STATVAR.vegetation.mlcanopyinst.qref; %forc_q(c); in kg/kg according to Bonan
    vegetation.mlcanopyinst.pref(p) = forcing.TEMP.p; % 101325.0; %  ground.STATVAR.vegetation.mlcanopyinst.pref; %forc_pbot(c);
    vegetation.mlcanopyinst.irsky(p) = forcing.TEMP.Lin; %ground.STATVAR.vegetation.mlcanopyinst.irsky; %forc_lwrad(c);
    vegetation.mlcanopyinst.qflx_rain(p) = forcing.TEMP.rainfall ./ (24*3600); %ground.STATVAR.vegetation.mlcanopyinst.qflx_rain; %forc_rain(c);
    vegetation.mlcanopyinst.qflx_snow(p) = forcing.TEMP.snowfall ./ (24*3600); %ground.STATVAR.vegetation.mlcanopyinst.qflx_snow; %forc_snow(c);
    
    % CO2 and O2: note unit conversion
    
    vegetation.mlcanopyinst.co2ref (p) = vegetation.mlcanopyinst.co2ref; %forc_pco2(g) / forc_pbot(c) * 1.e06;   % Pa -> umol/mol
    vegetation.mlcanopyinst.o2ref(p)  = vegetation.mlcanopyinst.o2ref; %forc_po2(g) / forc_pbot(c) * 1.e03;    % Pa -> mmol/mol
    
    % Miscellaneous
    %
    %     btran(p) = btran_patch(p);
    %     tacclim(p) = t_a10_patch(p);
    
    % Check to see if forcing height has changed
    
    if (vegetation.mlcanopyinst.zref(p) ~= vegetation.mlcanopyinst.zref_old(p))
        error (' ERROR: CanopyFluxesMultilayer: forcing height is not constant');
    end
    vegetation.mlcanopyinst.zref_old(p) = vegetation.mlcanopyinst.zref(p);
    
end

%---------------------------------------------------------------------
% Derived atmospheric input
%---------------------------------------------------------------------
for f = 1:vegetation.canopy.num_exposedvegp
    p = vegetation.canopy.filter_exposedvegp(f);
    vegetation.mlcanopyinst.eref(p) = vegetation.mlcanopyinst.qref(p) * vegetation.mlcanopyinst.pref (p) / (vegetation.physcon.mmh2o / vegetation.physcon.mmdry + (1  - vegetation.physcon.mmh2o / vegetation.physcon.mmdry) * vegetation.mlcanopyinst.qref(p));
    %     vegetation.mlcanopyinst.rhomol(p) = vegetation.mlcanopyinst.pref(p) / (vegetation.physcon.rgasc * vegetation.mlcanopyinst.tref(p));
    vegetation.mlcanopyinst.rhomol(p) = vegetation.mlcanopyinst.pref(p) / (8.3145 * vegetation.mlcanopyinst.tref(p));
    vegetation.mlcanopyinst.rhoair(p) = vegetation.mlcanopyinst.rhomol(p) * vegetation.physcon.mmdry * (1  - (1  - vegetation.physcon.mmh2o/vegetation.physcon.mmdry) * vegetation.mlcanopyinst.eref(p) / vegetation.mlcanopyinst.pref (p));
    vegetation.mlcanopyinst.mmair(p) = vegetation.mlcanopyinst.rhoair(p) / vegetation.mlcanopyinst.rhomol(p);
    vegetation.mlcanopyinst.cpair(p) = vegetation.physcon.cpd * (1  + (vegetation.physcon.cpw/vegetation.physcon.cpd - 1 ) * vegetation.mlcanopyinst.qref(p)) * vegetation.mlcanopyinst.mmair(p); %1200; %
    vegetation.mlcanopyinst.thref(p) = vegetation.mlcanopyinst.tref(p) + 0.0098  * vegetation.mlcanopyinst.zref(p);
    vegetation.mlcanopyinst.thvref(p) = vegetation.mlcanopyinst.thref(p) * (1  + 0.61  * vegetation.mlcanopyinst.qref(p));
end