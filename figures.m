function vegetation = figures(vegetation)


figure (1)
plot([vegetation.mlcanopyinst.tleaf(:,:,1)'-273.15 vegetation.mlcanopyinst.tair(:,:)'-273.15 vegetation.mlcanopyinst.tveg(:,:,1)'-273.15])
axis([1 25 -20 30])
% set(gcf,'position',[300 50 400 250])
% set(gcf,'position',[360 200 600 500])
title('Temperature of a diurnal cycle')
xlabel('Canopy layer')
ylabel('Temperature (°C)')
set(gcf,'position',[2840 800 600 500])

pause(0.1)
hold on

T = table(vegetation.mlcanopyinst.tleaf(:,:,1)'-273.15, vegetation.mlcanopyinst.tair(:,:)'-273.15, 'VariableNames', { 'A', 'B'} );
% Write data to text file
writetable(T, 'MyFile.txt')

% Create a table with the data and variable names


% plot(vegetation.mlcanopyinst.tg'-273.15,'*')



figure (2) %blue red yellow purple
plot(vegetation.mlcanopyinst.st_prof')
axis([1 25 -10 100])
% set(gcf,'position',[300 200 600 450])
% set(gcf,'position',[360 800 600 500])
title('Storage heat flux')
set(gcf,'position',[2220 200 600 500])
hold on
pause(0.1)
% plot(vegetation.mlcanopyinst.gsoi,'*');

figure (3)
plot(vegetation.mlcanopyinst.sh_prof')
axis([1 25 -40 100])
% set(gcf,'position',[980 200 600 500])
title('Sensibel heat fluxes of a diurnal cycle')
xlabel('Canopy layer')
ylabel('Flux (W/m^2)')
set(gcf,'position',[1600 200 600 500])
hold on
pause(0.1)

% canopy photosynthesis profile (last timestep)
figure (9)
plot(vegetation.mlcanopyinst.an(:,:,1))
plot(vegetation.mlcanopyinst.an(:,:,2))
set(gcf,'position',[2840 200 600 500])
title('Leaf photosynthesis rate profile') % https://faculty.washington.edu/edford/acclimation.html and http://www.jayreimer.com/TEXTBOOK/iText/products/0-13-115516-4/ch8/ch8_s3_5.html
hold on
ylabel('Photosynthesis rate (umol CO2/m^2 leaf/s)')
xlabel('Canopy layer')
pause(0.1)

figure (4)
plot(vegetation.mlcanopyinst.lh_prof')
axis([1 25 -40 100])
% set(gcf,'position',[980 800 600 500])
title('Latent heat fluxes of a diurnal cycle')
xlabel('Canopy layer')
ylabel('Flux (W/m^2)')
set(gcf,'position',[900 200 600 500])
hold on
pause(0.1)

figure (5)
plot (vegetation.mlcanopyinst.rn_prof')
axis([1 25 -40 150])
% set(gcf,'position',[1600 200 600 500])
title('Net radiation of a diurnal cycle')
xlabel('Canopy layer')
ylabel('Flux (W/m^2)')
set(gcf,'position',[200 200 600 500],'DefaultFigureVisible','off')
hold on
pause(0.1)
% plot(vegetation.flux.swsoi','*')
% plot(vegetation.flux.irsoi','x')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figures
height = vegetation.mlcanopyinst.zw./vegetation.mlcanopyinst.ztop;

figure (6)
% wind profile
plot(vegetation.mlcanopyinst.wind,height);
set(gcf,'position',[1600 800 600 500])
hold on
pause(0.1)
title('Wind profile (last timestep)')
xlabel('Windspeed (m/s^2)')
ylabel('Fraction of canopy height')

% figure (8)
% % canopy temperature profile (last timestep)
% plot(vegetation.mlcanopyinst.tair(:,:,1)-273.15,height)
% % plot(vegetation.mlcanopyinst.tleaf(:,:,1)-273.15,height)
% % plot(vegetation.mlcanopyinst.tleaf(:,:,2)-273.15,height)
% plot(vegetation.mlcanopyinst.tveg(:,:,2)-273.15,height)
% set(gcf,'position',[2220 800 600 500])
% title('Temperature profile (last timestep)')
% xlabel('Temperature (°C)')
% ylabel('Fraction of canopy height')


figure (11)
% canopy shaded / sunny fraction
plot(vegetation.flux.fracsun,height)
plot(vegetation.flux.fracsha,height)
set(gcf,'position',[900 800 600 500])
hold on
pause(0.1)
title('Fraction of sunny and shaded leaves')
xlabel('Fraction sunny (blue) vs. shaded (red)')
ylabel('Fraction of canopy height')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





% figure (10)
% plot(vegetation.mlcanopyinst.ustar', '*')
% title('Friction velocity u*')
% hold on
% % set(gcf,'position',[2200 800 600 500])
% pause(0.1)



end
