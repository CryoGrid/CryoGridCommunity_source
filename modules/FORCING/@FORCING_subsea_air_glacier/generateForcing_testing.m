function forcing = generateForcing_testing(forcing)
    timeForcing = [forcing.PARA.startForcing:forcing.PARA.dtForcing:forcing.PARA.endForcing]';
    forcing.DATA.timeForcing = timeForcing;
    forcingData = zeros(size(forcing.DATA.timeForcing));
    saltConcForcing  = zeros(size(forcing.DATA.timeForcing));
    coverage = zeros(size(forcing.DATA.timeForcing));
    time_inundation = 0;
    surfaceState = zeros(size(forcing.DATA.timeForcing));

    testcase = forcing.PARA.TF_no;

switch testcase
    case 1
        %one-third glacier, one third aerial, one third flooded

        i_start = 1;
        i_end = floor(length(forcingData)/3);
        forcingData(i_start:i_end) = forcing.PARA.T_IceSheet;
        coverage(i_start:i_end) = forcing.PARA.IS +100;
        surfaceState(i_start:i_end) = -1;

        i_start = i_end + 1;
        i_end = 2*i_end;
        forcingData(i_start:i_end) = linspace(-17, -9, i_end - i_start+1);
        coverage(i_start:i_end) = 0;
        surfaceState(i_start:i_end) = 1;

        i_start = i_end + 1;
        i_end = length(forcingData);
        forcingData(i_start:i_end) = linspace(forcing.PARA.T_freeze, 0, i_end - i_start+1);
        coverage(i_start:i_end) = linspace(-30, -2, i_end - i_start+1);
        surfaceState(i_start:i_end) = 0;
        saltConcForcing(i_start:i_end) = forcing.PARA.benthicSalt;

        time_inundation = i_end - i_start+1;

    case 2 %subaerial ramp
        forcingData = linspace(-20,20, length(forcing.DATA.timeForcing));
        coverage = zeros(size(forcing.DATA.timeForcing));
        time_inundation = 0;
        surfaceState = ones(size(forcing.DATA.timeForcing));

    case 3 %zero temperature, saltConc change with ramp

        i_start = floor(length(forcingData)/2);
        forcingData = 0.*ones(size(forcing.DATA.timeForcing));
        coverage = zeros(size(forcing.DATA.timeForcing));
        coverage(i_start:end) = 20;
        time_inundation = length(forcing.DATA.timeForcing) - i_start;
        surfaceState = ones(size(forcing.DATA.timeForcing));
        surfaceState(i_start:end) = 0;

        i_start = floor(length(forcingData)/3);
        i_end = 2*i_start;
        saltConcForcing(i_start:i_end) = linspace(0,forcing.PARA.benthicSalt, i_end - i_start + 1);
        saltConcForcing(i_end:end) = forcing.PARA.benthicSalt;

    case 4 %subsea ramp

        forcingData = linspace(5,-10, length(forcing.DATA.timeForcing));
        coverage = linspace(2,200, length(forcing.DATA.timeForcing));
        time_inundation = length(forcing.DATA.timeForcing);
        surfaceState = ones(size(forcing.DATA.timeForcing));

        saltConcForcing(:) = forcing.PARA.benthicSalt;
            
    case 5 %zero temperature, saltConc change without ramp

        i_start = floor(length(forcingData)/2);
        forcingData = -10.*ones(size(forcing.DATA.timeForcing));
        coverage = zeros(size(forcing.DATA.timeForcing));
        coverage(i_start:end) = 20;
        time_inundation = length(forcing.DATA.timeForcing) - i_start;
        surfaceState = ones(size(forcing.DATA.timeForcing));
        surfaceState(i_start:end) = 0;

        saltConcForcing(i_start:end) = forcing.PARA.benthicSalt;
        
            
    case 6 %rapid temperature and saltConc change

        i_start = floor(length(forcingData)/6);
        
        saltConcForcing(1:i_start) = 0;
        saltConcForcing(i_start:end) = forcing.PARA.benthicSalt;

        i_start = 2*i_start;
        forcingData(1:i_start) = -10;
        forcingData(i_start:end) = 0;
end



forcing.DATA.TForcing = forcingData;
forcing.DATA.coverage = coverage;
forcing.DATA.time_inundation = time_inundation;
forcing.DATA.surfaceState = surfaceState;
if forcing.PARA.saltForcingSwitch == 0
    forcing.DATA.saltConcForcing = zeros(size(saltConcForcing));
else
    forcing.DATA.saltConcForcing = saltConcForcing;
end


% figure(4)
% plot(timeForcing, forcingData)
% hold on
% plot(timeForcing, saltConcForcing)
end
