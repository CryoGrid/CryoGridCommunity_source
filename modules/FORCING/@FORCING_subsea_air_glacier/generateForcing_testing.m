function forcing = generateForcing_testing(forcing)
    timeForcing = [forcing.PARA.startForcing:forcing.PARA.dtForcing:forcing.PARA.endForcing]';
    forcing.DATA.timeForcing = timeForcing;
    
    airTemp = zeros(size(forcing.DATA.timeForcing));
    seaLevel = zeros(size(forcing.DATA.timeForcing));
    glacialCover = zeros(size(forcing.DATA.timeForcing));


    testcase = forcing.PARA.TF_no;

switch testcase
    case 1
        %one-third glacier
        i_start = 1;
        i_end = floor(length(timeForcing)/3);
        glacialCover(i_start:i_end) = forcing.PARA.IS +100;

        %one third aerial
        i_start = i_end + 1;
        i_end = 2*i_end;
        airTemp(i_start:i_end) = linspace(-17, -9, i_end - i_start+1);

        %one third flooded
        i_start = i_end + 1;
        i_end = length(timeForcing);
        seaLevel(i_start:i_end) = linspace(30, 100, i_end - i_start+1);

    case 2 %subaerial ramp
        airTemp = linspace(-20,20, length(forcing.DATA.timeForcing));
        glacialCover = zeros(size(forcing.DATA.timeForcing));
        seaLevel = -100*ones(size(forcing.DATA.timeForcing));

    case 3 %zero temperature, saltConc change with ramp

        i_start = floor(length(timeForcing)/2);
        airTemp = 0.*ones(size(forcing.DATA.timeForcing));
        
        seaLevel(i_start:end) = 20;

    case 4 %subsea ramp

        airTemp = linspace(5,-10, length(forcing.DATA.timeForcing));
        seaLevel = linspace(5,100, length(forcing.DATA.timeForcing));
            
    case 5 %zero temperature, saltConc change without ramp

        i_start = floor(length(timeForcing)/2);
        airTemp = -10.*ones(size(forcing.DATA.timeForcing));
        
        seaLevel(i_start:end) = 20;
            
    case 6 %rapid temperature and saltConc change

        i_start = floor(length(timeForcing)/6);
        
        seaLevel(1:i_start) = -100;
        seaLevel(i_start:end) = 100;

        i_start = 2*i_start;
        airTemp(1:i_start) = -10;
        airTemp(i_start:end) = 0;
end



forcing.DATA.airTemp = airTemp;
forcing.DATA.seaLevel = seaLevel;
forcing.DATA.glacialCover = glacialCover;
end
