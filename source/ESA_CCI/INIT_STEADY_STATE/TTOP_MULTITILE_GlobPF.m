classdef TTOP_MULTITILE_GlobPF < BASE
    
    %own class that can be called by TILE, using and rearranging
    %FORCING.DATA, or reading data from elsewwhere
    

    methods
        
        function ground = provide_PARA(ground)
            ground.PARA.start_year = [];
            ground.PARA.end_year = [];
            ground.PARA.nf_folder = [];
            
            ground.PARA.rk_init = []; %from ensemble
            ground.PARA.snowfall_factor = []; 
            ground.PARA.bare_forest_fraction = [];
        end
        
        function ground = provide_CONST(ground)
            
        end
        
        function ground = provide_STATVAR(ground)
            
        end
        
        function ground = finalize_init(ground, tile)
            
        end
        
        function TTOP = get_TTOP_from_forcing(ground, tile)
            
            %(FORCING, year_list, snowfall_factor, ensemble_size, rk) %
            timestamp = tile.FORCING.DATA.timestamp;
            start_pos=1;
            while str2num(datestr(timestamp(1,start_pos), 'yyyy')) < ground.PARA.start_year
                start_pos=start_pos+1;
            end
            end_pos=1;
            while str2num(datestr(timestamp(1,end_pos), 'yyyy')) <= ground.PARA.end_year
                end_pos=end_pos+1;
            end
            end_pos = end_pos-1;
            shift_pos=start_pos;
            while str2num(datestr(timestamp(1,shift_pos), 'mm'))< 8
                shift_pos=shift_pos+1;
            end
            shift_yes = tile.PARA.latitude > 0; %Northern hemisphere
            
            timestamp = timestamp(:, start_pos:end_pos);
            surf_T = tile.FORCING.DATA.final_av_T(:, start_pos:end_pos);
            snowfall = tile.FORCING.DATA.ERA_snowfall_downscaled(:, start_pos:end_pos);
            melt_bare = tile.FORCING.DATA.ERA_melt_bare(:, start_pos:end_pos);
            melt_forest = tile.FORCING.DATA.ERA_melt_forest(:, start_pos:end_pos);

            January = find(str2num(datestr(timestamp, 'mm'))==1);
            July = find(str2num(datestr(timestamp, 'mm'))==7); %Southern Hemisphere
            Jan_av_T = mean(surf_T(:,January),2);
            Jan_av_T(~shift_yes,1) = mean(surf_T(~shift_yes, July),2);
            
            surf_T_shift = [tile.FORCING.DATA.final_av_T(shift_yes,shift_pos:end_pos) tile.FORCING.DATA.final_av_T(shift_yes,start_pos:shift_pos-1)];
            snowfall_shift = [tile.FORCING.DATA.ERA_snowfall_downscaled(shift_yes,shift_pos:end_pos) tile.FORCING.DATA.ERA_snowfall_downscaled(shift_yes,start_pos:shift_pos-1)];
            melt_bare_shift = [tile.FORCING.DATA.ERA_melt_bare(shift_yes,shift_pos:end_pos) tile.FORCING.DATA.ERA_melt_bare(shift_yes,start_pos:shift_pos-1)];
            melt_forest_shift = [tile.FORCING.DATA.ERA_melt_forest(shift_yes,shift_pos:end_pos) tile.FORCING.DATA.ERA_melt_forest(shift_yes,start_pos:shift_pos-1)];
            %timestamp_shift = [tile.FORCING.DATA.timestamp(1,shift_pos:end_pos) tile.FORCING.DATA.timestamp(1,start_pos:shift_pos-1)];
            
            surf_T(shift_yes,:) = surf_T_shift;
            snowfall(shift_yes,:) = snowfall_shift;
            melt_bare(shift_yes,:) = melt_bare_shift;
            melt_forest(shift_yes,:) = melt_forest_shift;
            
%             FDD = repmat(FDD', 1, tile.ENSEMBLE.PARA.ensemble_size);
%             TDD = repmat(TDD', 1, tile.ENSEMBLE.PARA.ensemble_size);
            Jan_av_T = repmat(Jan_av_T', 1, tile.PARA.ensemble_size);
            
            snow_cover_duration_av = 0;
                        
            FDD_acc=0;
            TDD_acc = 0;
            swe=0;
            swe_acc=0;
            swe_number=0;
            swe_test=[];
            for i=1:size(snowfall,2)
                T = repmat(surf_T(:,i)', 1, tile.PARA.ensemble_size);
                sf = repmat(snowfall(:,i)', 1, tile.PARA.ensemble_size) .* ground.PARA.snowfall_factor;
                melt = repmat(melt_bare(:,i)', 1, tile.PARA.ensemble_size) .* (1 - ground.PARA.bare_forest_fraction) + ...
                    repmat(melt_forest(:,i)', 1, tile.PARA.ensemble_size) .* ground.PARA.bare_forest_fraction;
                
                swe = swe + (sf - melt) .* (timestamp(1,2) - timestamp(1,1));
                swe(swe<0)=0;
                swe_acc = swe_acc+swe;
                swe_number = swe_number + double(swe>0);
                snow_cover_duration_av = snow_cover_duration_av + double(swe>0);
                FDD_acc = FDD_acc + double(T<0) .* T;
                TDD_acc = TDD_acc + double(T>0 & swe<=0) .* T;
            end
            
            swe_av = swe_acc ./ swe_number;
            
            snow_cover_duration_av = snow_cover_duration_av ./ size(snowfall,2) .* 12;  %CHANGED
            
            TTOP = get_GlobPF_TTOP(ground, FDD_acc, TDD_acc, swe_av ./ 1000, Jan_av_T, snow_cover_duration_av, ground.PARA.rk_init, size(snowfall,2));
            
            %REMOVE THIS!
            ground.STATVAR.snow_cover_duration_av = snow_cover_duration_av;
            ground.STATVAR.swe_av = swe_av;
            ground.STATVAR.Jan_av_T = Jan_av_T;
            ground.STATVAR.FDD = FDD_acc;
            ground.STATVAR.TDD = TDD_acc;
        end
        
        function TTOP = get_GlobPF_TTOP(ground, FDD, TDD, av_SWE, MAAT_Jan, snow_cover_duration, rk, number_of_days)
            
            load([ground.PARA.nf_folder 'nf4GlobPF.mat']);
            
            
            MAST=(FDD+TDD)./number_of_days;  %mean annual surface temperature
            
            %--------------------------
            duration_snow_cover=max(snow_cover_duration,1);
            MAAT_Jan=min(MAAT_Jan, -2); %density values become a bit unphysical here
            
            %based on Onuchin and Burenia (1996)
            snow_density=1000.*(0.5.*(0.323-0.0477.*log(-MAAT_Jan)) + sqrt( (0.5.*(0.323-0.0477.*log(-MAAT_Jan))).^2 + 0.00045.*0.917.*av_SWE./10.*log(duration_snow_cover))); % in kg/m3
            
            
            snow_depth=av_SWE.*917./snow_density;
            
            MAST=max(MAST.*0+min(MAAT), min(MAST.*0+max(MAAT), MAST));
            snow_depth=max(snow_depth.*0+min(sd), min(snow_depth.*0+max(sd), snow_depth));
            
            
            nf=interp2(MAAT, sd, nf_param, MAST, snow_depth);
            
            TTOP = double(nf.*FDD + rk.*TDD <= 0) .* (nf.*FDD + rk.*TDD)./number_of_days + double(nf.*FDD + rk.*TDD > 0) .* (nf./rk.*FDD + TDD)./number_of_days;
            
            
        end
        
        function TTOP = get_TTOP(ground, FDD, TDD, rk, number_of_days)
            
            TTOP = double(FDD + rk.*TDD <= 0) .* (FDD + rk.*TDD)./number_of_days + double(FDD + rk.*TDD > 0) .* (1./rk.*FDD + TDD)./number_of_days;
            
        end
        
    end
end