function ok = interpolateResultsSubseaPF( input_args )
%INTERPOLATERESULTSSUBSEAPF Summary of this function goes here
%   Detailed explanation goes here

% T=zeros(length(cT_grid),length(modelTimeSlice));  % skip ODE call...

%% replace grid cells above the upper boundary by NaN
for i=1:size(modelTimeSlice,1)
   ub=interp1(uppermostGridCell(:,1), uppermostGridCell(:,2), modelTimeSlice(i), 'nearest');
   if ub>1
        T(1:ub-1,i)=NaN;
   end
end
 
%% output diagnostics
pars=zeros(1000,10);

%% interpolate temperature & parameters from 6 km-cT_grid to 2 km-2 m grid
TS = interp2(modelTimeSlice, cT_grid, T, modelTimeSlice, [1:2:1999]); % soil temperature on grid
% pars:   1: cT_grid;  2: mineral;  3: water;   4: organic;   
%         5: salinity;   6: soilType;   7: a;   8: b;  9: T_melt;   10: terrestrialYesNo;

%Interpolate Linear:
%1: cT_grid;  2: mineral;  3: water;   4: organic;   
for ii = 1:4
    pars(:,ii) = interp1(cT_grid, parameters(:,ii),[1:2:1999],'linear'); 
end

%Interpolate with nearest neighbour:
%5: salinity;   6: soilType;   7: a;   8: b;  9: T_melt;   10: terrestrialYesNo;
for ii = 5:10
    pars(:,ii) = interp1(cT_grid, parameters(:,ii),[1:2:1999],'neareast');
end

%% calculate ice content
waterSaturation=[];
iceSaturation=[];
liquidWC=[];
for i=1:size(modelTimeSlice,1)
    %a1=1./pars(:,3) - pars(:,7).*pars(:,9) - pars(:,8).*pars(:,9).^2;
    %licWC = (TS(:,i)<pars(:,9)).*(1./(a1 + pars(:,7).*TS(:,i) + pars(:,8).*TS(:,i).^2)) + (TS(:,i)>=pars(:,9)).*pars(:,3);
    licWC = get_wc(T(:,i), data);
    licWC = interp1(cT_grid,licWC, [1:2:1999],'linear')';
    liquidWC = [liquidWC licWC];
    iceSaturation=[iceSaturation (pars(:,3)-licWC)./pars(:,3)];
    waterSaturation=[waterSaturation licWC./pars(:,3)];
end

%% find upper and lower PF boundaries based on soil temperatures TS
PF_UB = zeros(length(modelTimeSlice),1);  % permafrost lower boundary
PF_LB = ones(length(modelTimeSlice),1)*1999;  % permafrost lower boundary
for j = 1:length(modelTimeSlice) % time loop
    i=1; % upper PF boundary
    while TS(i,j)>0  && i<size(TS,1) % z loop
        PF_UB(j)=PF_UB(j)+2;
        i=i+1;
    end
    PF_UB(j)=PF_UB(j)-1;

    i=size(TS,1); % lower PF boundary
    while TS(i,j)>0  && i>1
        PF_LB(j)=PF_LB(j)-2;
        i=i-1;
    end
    PFB = [PF_UB PF_LB];
end

%% write output files to structure variable

resstruc = struct('time',single(timestamp(2:end-1)), ...
                 'modelTime', single(modelTimeSlice), ...
                 'timeF',single(Tsurf(2:end-1,1)), ...
                 'TF',single(Tsurf(2:end-1,2)), ...
                 'TS',single(TS), 'PFB',single(PFB), ...
                 'IceSat',single(iceSaturation), ...
                 'surfState',single(surfaceState), ...
                 'coverage', single(coverage), ...
                 'Lat', single(lat), 'Lon', single(lon), ...
                 'pars',single(pars), ...
                 'EaseGridNr', single(nr), 'EaseGridRes', single(EDF.GR)); 

        
%save results with their EaseGrid Number or LatLon Values
filename=strcat(outputFolderName,filesep, sprintf('resultstruc_Lat%03.6f_Lon%03.6f', lat, lon), '.mat');
try
    save(filename,'resstruc') 
end

