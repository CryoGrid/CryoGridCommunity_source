% figures after parallel runs

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 1
% % Ground under vegetation (winter)
% 
% % Ground depth: 
% depth = cumsum(out.STRATIGRAPHY{1, 1}{2, 1}.STATVAR.layerThick);
% Temp = [];
% Temp_snow = [];
% for i = 1:80
% %     [numRows,numCols] = size(out.STRATIGRAPHY{1, i});
% %     if numRows == 2
%     out.STRATIGRAPHY{1, i}{2, 1}.STATVAR.T(:,2) = depth(:);
%     Temp(:,i) = out.STRATIGRAPHY{1, i}{2, 1}.STATVAR.T(:,1);
% %     else
% %     Temp_snow(:,i) = out.STRATIGRAPHY{1, i}{3, 1}.STATVAR.T(:,1);
% %     end    
% end
% 
% nanmean(Temp(1,:));
% 
% for o = 1:450
% mean_t(o,:) = nanmean(Temp(o,:));
% median_t(o,:) = nanmedian(Temp(o,:));
% min_t(o,:) = nanmin(Temp(o,:));
% max_t(o,:) = nanmax(Temp(o,:));
% % mean_t_snow(o,:) = nanmean(Temp_snow(o,:));
% % median_t_snow(o,:) = nanmedian(Temp_snow(o,:));
% % max_t_snow(o,:) = nanmax(Temp_snow(o,:));
% % min_t_snow(o,:) = nanmin(Temp_snow(o,:));
% end
% 
% figure (1) 
% plot( mean_t(1:200,1), depth(1:200,1));
% hold on
% plot(median_t(1:200,1), depth(1:200,1));
% plot(min_t(1:200,1), depth(1:200,1));
% plot(max_t(1:200,1), depth(1:200,1));
% % plot(depth(1:200,1), mean_t_snow(1:200,1),':');
% % plot(depth(1:200,1), median_t_snow(1:200,1),':');
% % plot(depth(1:200,1), min_t_snow(1:200,1),':');
% % plot(depth(1:200,1), max_t_snow(1:200,1),':');
% title('Ground temperature')
% ylabel('Ground depth (1-12 m)')
% axis ij
% xlabel('Temperature (°C)')
% axis([-20 20 0 12])
% 
% 
% % Vegetation temperature: 
% height_vegetation = out.STRATIGRAPHY{1, 1}{1, 1}.STATVAR.vegetation.mlcanopyinst.zw./out.STRATIGRAPHY{1, 1}{1, 1}.STATVAR.vegetation.mlcanopyinst.ztop;
% Temp_vegetation = [];
% for i = 1:364
%     [numRows,numCols] = size(out.STRATIGRAPHY{1, i});
%     if numRows == 2
%     Temp_vegetation(:,i) = out.STRATIGRAPHY{1, i}{1, 1}.STATVAR.vegetation.mlcanopyinst.tveg(:,:,1)-273.15;
%     else
%     Temp_vegetation_snow(:,i) = out.STRATIGRAPHY{1, i}{2, 1}.STATVAR.vegetation.mlcanopyinst.tveg(:,:,1)-273.15;
%     end    
% end
% 
% nanmean(Temp_vegetation(1,:));
% 
% for o = 1:15
% mean_t_vegetation(o,:) = nanmean(Temp_vegetation(o,:));
% median_t_vegetation(o,:) = nanmedian(Temp_vegetation(o,:));
% min_t_vegetation(o,:) = nanmin(Temp_vegetation(o,:));
% max_t_vegetation(o,:) = nanmax(Temp_vegetation(o,:));
% mean_t_vegetation_snow(o,:) = nanmean(Temp_vegetation_snow(o,:));
% median_t_vegetation_snow(o,:) = nanmedian(Temp_vegetation_snow(o,:));
% max_t_vegetation_snow(o,:) = nanmax(Temp_vegetation_snow(o,:));
% min_t_vegetation_snow(o,:) = nanmin(Temp_vegetation_snow(o,:));
% end
% 
% figure (2) 
% plot(mean_t_vegetation(1:15,1), depth(1:15,1));
% hold on
% plot(median_t_vegetation(1:15,1), depth(1:15,1));
% plot(min_t_vegetation(1:15,1), depth(1:15,1));
% plot(max_t_vegetation(1:15,1), depth(1:15,1));
% plot(mean_t_vegetation_snow(1:15,1), depth(1:15,1), ':');
% plot(median_t_vegetation_snow(1:15,1), depth(1:15,1),':');
% plot(min_t_vegetation_snow(1:15,1), depth(1:15,1),':');
% plot(max_t_vegetation_snow(1:15,1), depth(1:15,1),':');
% title('Vegetation temperature')
% ylabel('Vegetation height (1-8 m)')
% xlabel('Temperature (°C)')
% axis([-20 20 0 1])

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 2
% % Ground under vegetation (summer)
% 
% % Ground depth: 
% depth = cumsum(out.STRATIGRAPHY{1, 1}{2, 1}.STATVAR.layerThick);
% Temp = [];
% % Temp_snow = [];
% for i = 1:80
% %     [numRows,numCols] = size(out.STRATIGRAPHY{1, i});
% %     if numRows == 2
%     out.STRATIGRAPHY{1, i}{2, 1}.STATVAR.T(:,2) = depth(:);
%     Temp(:,i) = out.STRATIGRAPHY{1, i}{2, 1}.STATVAR.T(:,1);
% %     else
% %     Temp_snow(:,i) = out.STRATIGRAPHY{1, i}{3, 1}.STATVAR.T(:,1);
% %     end    
% end
% 
% nanmean(Temp(1,:));
% 
% for o = 1:470
% mean_t(o,:) = nanmean(Temp(o,:));
% median_t(o,:) = nanmedian(Temp(o,:));
% min_t(o,:) = nanmin(Temp(o,:));
% max_t(o,:) = nanmax(Temp(o,:));
% % mean_t_snow(o,:) = nanmean(Temp_snow(o,:));
% % median_t_snow(o,:) = nanmedian(Temp_snow(o,:));
% % max_t_snow(o,:) = nanmax(Temp_snow(o,:));
% % min_t_snow(o,:) = nanmin(Temp_snow(o,:));
% end
% 
% figure (1) 
% plot( mean_t(1:200,1), depth(1:200,1));
% hold on
% plot(median_t(1:200,1), depth(1:200,1));
% plot(min_t(1:200,1), depth(1:200,1));
% plot(max_t(1:200,1), depth(1:200,1));
% % plot(depth(1:200,1), mean_t_snow(1:200,1),':');
% % plot(depth(1:200,1), median_t_snow(1:200,1),':');
% % plot(depth(1:200,1), min_t_snow(1:200,1),':');
% % plot(depth(1:200,1), max_t_snow(1:200,1),':');
% title('Ground temperature')
% ylabel('Ground depth (1-12 m)')
% axis ij
% xlabel('Temperature (°C)')
% axis([-10 40 0 12])
% 
% 
% % Vegetation temperature: 
% height_vegetation = out.STRATIGRAPHY{1, 1}{1, 1}.STATVAR.vegetation.mlcanopyinst.zw./out.STRATIGRAPHY{1, 1}{1, 1}.STATVAR.vegetation.mlcanopyinst.ztop;
% Temp_vegetation = [];
% for i = 1:372
%     [numRows,numCols] = size(out.STRATIGRAPHY{1, i});
% %     if numRows == 2
%     Temp_vegetation(:,i) = out.STRATIGRAPHY{1, i}{1, 1}.STATVAR.vegetation.mlcanopyinst.tveg(:,:,1)-273.15;
% %     else
% %     Temp_vegetation_snow(:,i) = out.STRATIGRAPHY{1, i}{2, 1}.STATVAR.vegetation.mlcanopyinst.tveg(:,:,1)-273.15;
% %     end    
% end
% 
% nanmean(Temp_vegetation(1,:));
% 
% for o = 1:15
% mean_t_vegetation(o,:) = nanmean(Temp_vegetation(o,:));
% median_t_vegetation(o,:) = nanmedian(Temp_vegetation(o,:));
% min_t_vegetation(o,:) = nanmin(Temp_vegetation(o,:));
% max_t_vegetation(o,:) = nanmax(Temp_vegetation(o,:));
% % mean_t_vegetation_snow(o,:) = nanmean(Temp_vegetation_snow(o,:));
% % median_t_vegetation_snow(o,:) = nanmedian(Temp_vegetation_snow(o,:));
% % max_t_vegetation_snow(o,:) = nanmax(Temp_vegetation_snow(o,:));
% % min_t_vegetation_snow(o,:) = nanmin(Temp_vegetation_snow(o,:));
% end
% 
% figure (2) 
% plot(mean_t_vegetation(1:15,1), depth(1:15,1));
% hold on
% plot(median_t_vegetation(1:15,1), depth(1:15,1));
% plot(min_t_vegetation(1:15,1), depth(1:15,1));
% plot(max_t_vegetation(1:15,1), depth(1:15,1));
% % plot(mean_t_vegetation_snow(1:15,1), depth(1:15,1), ':');
% % plot(median_t_vegetation_snow(1:15,1), depth(1:15,1),':');
% % plot(min_t_vegetation_snow(1:15,1), depth(1:15,1),':');
% % plot(max_t_vegetation_snow(1:15,1), depth(1:15,1),':');
% title('Vegetation temperature')
% ylabel('Vegetation height (1-8 m)')
% xlabel('Temperature (°C)')
% axis([-10 40 0 1])

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3
% Ground without vegetation (summer)

% % Ground depth: 
% depth = cumsum(out.STRATIGRAPHY{1, 1}{1, 1}.STATVAR.layerThick);
% Temp = [];
% % Temp_snow = [];
% for i = 1:80
% %     [numRows,numCols] = size(out.STRATIGRAPHY{1, i});
% %     if numRows == 2
%     out.STRATIGRAPHY{1, i}{1, 1}.STATVAR.T(:,2) = depth(:);
%     Temp(:,i) = out.STRATIGRAPHY{1, i}{1, 1}.STATVAR.T(:,1);
% %     else
% %     Temp_snow(:,i) = out.STRATIGRAPHY{1, i}{3, 1}.STATVAR.T(:,1);
% %     end    
% end
% 
% nanmean(Temp(1,:));
% 
% for o = 1:470
% mean_t(o,:) = nanmean(Temp(o,:));
% median_t(o,:) = nanmedian(Temp(o,:));
% min_t(o,:) = nanmin(Temp(o,:));
% max_t(o,:) = nanmax(Temp(o,:));
% end
% 
% figure (1) 
% plot( mean_t(1:200,1), depth(1:200,1));
% hold on
% plot(median_t(1:200,1), depth(1:200,1));
% plot(min_t(1:200,1), depth(1:200,1));
% plot(max_t(1:200,1), depth(1:200,1));
% title('Ground temperature')
% ylabel('Ground depth (1-12 m)')
% axis ij
% xlabel('Temperature (°C)')
% axis([-10 40 0 12])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1
% Ground (winter)

% Ground depth: 
depth = cumsum(out.STRATIGRAPHY{1, 1}{1, 1}.STATVAR.layerThick);
Temp = [];
Temp_snow = [];
for i = 1:1460
    [numRows,numCols] = size(out.STRATIGRAPHY{1, i});
    if numRows == 1
        out.STRATIGRAPHY{1, i}{1, 1}.STATVAR.T(:,2) = depth;
        Temp(:,i) = out.STRATIGRAPHY{1, i}{1, 1}.STATVAR.T(:,1);
    else
        out.STRATIGRAPHY{1, i}{2, 1}.STATVAR.T(:,2) = depth;
        Temp_snow(:,i) = out.STRATIGRAPHY{1, i}{2, 1}.STATVAR.T(:,1);
        %     Temp_snow(:,i) = out.STRATIGRAPHY{1, i}{3, 1}.STATVAR.T(:,1);
    end
end

for o = 1:200
mean_t(o,:) = nanmean(Temp(o,:));
median_t(o,:) = nanmedian(Temp(o,:));
min_t(o,:) = nanmin(Temp(o,:));
max_t(o,:) = nanmax(Temp(o,:));
mean_t_snow(o,:) = nanmean(Temp_snow(o,:));
median_t_snow(o,:) = nanmedian(Temp_snow(o,:));
max_t_snow(o,:) = nanmax(Temp_snow(o,:));
min_t_snow(o,:) = nanmin(Temp_snow(o,:));
end

figure (1) 
plot( mean_t(1:200,1), depth(1:200,1));
hold on
plot(median_t(1:200,1), depth(1:200,1));
plot(min_t(1:200,1), depth(1:200,1));
plot(max_t(1:200,1), depth(1:200,1));
plot(depth(1:200,1), mean_t_snow(1:200,1),':');
plot(depth(1:200,1), median_t_snow(1:200,1),':');
plot(depth(1:200,1), min_t_snow(1:200,1),':');
plot(depth(1:200,1), max_t_snow(1:200,1),':');
title('Ground temperature')
ylabel('Ground depth (1-12 m)')
axis ij
xlabel('Temperature (°C)')
axis([-20 20 0 12])






