function display_ground_SubseaPFwSalt(filepath, plotdepth)
% to read and display temperatures from file

%in the moment, this interpolates the results to a "master grid"
%think about a better was to display results - maybe surf or patch?

if nargin <1
    %filepath = 'Results/test3_initial_stratigraphy/Results_test3_initial_stratigraphy.mat';
    %filepath = 'Results/test1_onlyMarine/Results_test1_onlyMarine.mat';
    %filepath = 'Results/test4_onlyTerr_wSalt/Results_test4_onlyTerr_wSalt.mat';
    %filepath = 'results/submarine_Benchmark_Testlocation1/Results_submarine_Benchmark_Testlocation1.mat';
    %filepath = 'results/submarine_PF_test_salt/Results_submarine_PF_test_salt.mat';
    filepath = 'results/submarine_PF_test_salt_all450k/Results_submarine_PF_test_salt_all450k.mat';
end
if nargin < 2
    plotdepth = -500;
end
addpath(genpath('modules'))
out = OUT_subseaPF;
load(filepath)


Time = out.TIMESTAMP./365.25;
master_midpoints = -[1:2:1999]; %think about this. Needs to adapt somehow to the upperPos of the first module

TForcing = out.FORCING.TForcing;
saltConcForcing = out.FORCING.saltConcForcing;
timeForcing = out.FORCING.timeForcing;

T = zeros(length(master_midpoints), length(Time));
saltConc = zeros(length(master_midpoints), length(Time));
thermCond = zeros(length(master_midpoints), length(Time));
c_eff = zeros(length(master_midpoints), length(Time));
liqWater = zeros(length(master_midpoints), length(Time));
totalSalt = zeros(1, length(Time));

for i=1:length(Time) %loop over saved time steps
    %Read out current temperature and grid
    RES_T=[];
    RES_saltConc=[];
    RES_layerThick=[];
    RES_thermCond = [];
    RES_c_eff = [];
    RES_liqWater=[];

    for j=size(out.STRATIGRAPHY{1,i},1):-1:1 %loop over all cells in the stratigraphy
        RES_T= [out.STRATIGRAPHY{1,i}{j,1}.STATVAR.T; RES_T];
        RES_saltConc= [out.STRATIGRAPHY{1,i}{j,1}.STATVAR.saltConc; RES_saltConc];
        RES_liqWater= [out.STRATIGRAPHY{1,i}{j,1}.STATVAR.liqWater; RES_liqWater];
        RES_layerThick= [out.STRATIGRAPHY{1,i}{j,1}.STATVAR.layerThick; RES_layerThick];
        RES_thermCond= [out.STRATIGRAPHY{1,i}{j,1}.STATVAR.thermCond; RES_thermCond];
        RES_c_eff= [out.STRATIGRAPHY{1,i}{j,1}.STATVAR.c_eff; RES_c_eff];        
    end
    
    RES_layerDepth = [out.STRATIGRAPHY{1,i}{j,1}.STATVAR.upperPos; out.STRATIGRAPHY{1,i}{j,1}.STATVAR.upperPos - cumsum(RES_layerThick)];
    RES_midptDepth = out.STRATIGRAPHY{1,i}{j,1}.STATVAR.upperPos - RES_layerThick(1)/2 - cumsum(RES_layerThick);
    
    totalSalt(i) = sum(RES_saltConc.*RES_liqWater.*RES_layerThick);
%     RES_grid = [0; cumsum(RES_grid)];
%     RES_grid =  RES_grid(end)-RES_grid;
%     RES_grid = [RES_grid(1); 0.5.*(RES_grid(2:end) + RES_grid(1:end-1))] + out.STRATIGRAPHY{1,i}{end,1}.STATVAR.lowerPos;
    T(:,i) = interp1(RES_midptDepth, RES_T, master_midpoints);
    saltConc(:,i) = interp1(RES_midptDepth, RES_saltConc, master_midpoints);
    liqWater(:,i) = interp1(RES_midptDepth, RES_liqWater, master_midpoints);
    [RES_layerDepth, I] = unique(RES_layerDepth);
    RES_thermCond = RES_thermCond(I);
    thermCond(:,i) = interp1(RES_layerDepth, RES_thermCond, master_midpoints);
    c_eff(:,i) = interp1(RES_midptDepth, RES_c_eff, master_midpoints);
end
figure('Position', [1 1 1810 1227])
clf

subplot(7,1,2)
surf(Time, master_midpoints, T, 'EdgeColor', 'none')
view(2)
hold on
contour(Time, master_midpoints, T, [0,0], 'color', 'black')
colorbar
axis xy
xlim([timeForcing(1), timeForcing(end)])
ylim([plotdepth,0])
title('Temperature over time')
ylabel('depth / m')

subplot(7,1,3)
surf(Time, master_midpoints, saltConc, 'EdgeColor', 'none')
view(2)
ylim([plotdepth,0])
set(gca, 'YDir', 'reverse')
colorbar
axis xy
xlim([timeForcing(1), timeForcing(end)])
title('Salt Concentration over time')
ylabel('depth / m')

%get axes position with colorbar
pos0 = get(gca, 'Position');

subplot(7,1,4)
plot(Time, totalSalt);
xlim([timeForcing(1), timeForcing(end)])
pos = get(gca, 'Position');
set(gca, 'Position', [pos(1) pos(2) pos0(3) pos(4)]);
title('Total Salt over time')

subplot(7,1,5)
surf(Time, master_midpoints, liqWater, 'EdgeColor', 'none')
view(2)
ylim([plotdepth,0])
xlim([timeForcing(1), timeForcing(end)])
set(gca, 'YDir', 'reverse')
colorbar
axis xy
title('Liquid Water Content over time')
ylabel('depth / m')

subplot(7,1,6)
surf(Time, master_midpoints, thermCond, 'EdgeColor', 'none')
view(2)
xlim([timeForcing(1), timeForcing(end)])
ylim([plotdepth,0])
title('Thermal conductivity over time')
colorbar
axis xy
ylabel('depth / m')

subplot(7,1,7)
surf(Time, master_midpoints, c_eff, 'EdgeColor', 'none')
view(2)
xlim([timeForcing(1), timeForcing(end)])
ylim([plotdepth,0])
axis xy
colorbar
caxis([0,4e7])
title('Heat capacity over time')
ylabel('depth / m')
xlabel('time / a')

ax1 = subplot(7,1,1);
%yyaxis(ax1,'left')
ax1 = plotyy(timeForcing, TForcing, timeForcing, saltConcForcing);
ylabel(ax1(1), 'temperature / ^\circ C')
%yyaxis(ax1,'right')
hold on
%plot(timeForcing, saltConcForcing)
ylabel(ax1(2), 'salt concentration / mol / m^3')
title('Forcing Data over time')
xlim([timeForcing(1), timeForcing(end)])
pos = get(gca, 'Position');
set(gca, 'Position', [pos(1) pos(2) pos0(3) pos(4)]);


[path,name,ext] = fileparts(filepath);
saveas(gcf, fullfile(path, 'Results_overview.png'))


fprintf('Timesteps: %.0f \n', out.RUNINFO.timesteps);
fprintf('Minimal Timestep: %3.2f days \n', out.RUNINFO.dt_min);
fprintf('Maximum Timestep: %3.2f days \n', out.RUNINFO.dt_max);
fprintf('Runtime: %.2f hours \n', out.RUNINFO.endtime./(60*60));


% figure
% plot(TForcing)
% figure
% plot(D_surf')
% figure
% plot(density)

save(fullfile(path, 'Results_Matrices.mat'), 'Time', 'master_midpoints', 'TForcing', 'saltConcForcing', ...
                                              'timeForcing', 'T', 'saltConc', 'thermCond', 'c_eff', ...
                                              'liqWater', 'totalSalt');