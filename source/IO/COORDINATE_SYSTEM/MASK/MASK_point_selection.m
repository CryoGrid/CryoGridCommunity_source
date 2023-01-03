%========================================================================
% CryoGrid MASK class MASK_point_selection
% selects the region of interest as the closest coordinates to a 
% list of target coordinates provided in a file 
%
% S. Westermann, Jan 2021
% S. Westermann, Dec 2022
%========================================================================

classdef MASK_point_selection < matlab.mixin.Copyable

    
    properties
        PARENT
        PARA
        CONST
        STATVAR
    end
    
    methods
        
        function mask = provide_PARA(mask)
            mask.PARA.point_file_path = [];
            mask.PARA.filename = [];
            mask.PARA.additive = [];
        end

        function mask = provide_STATVAR(mask)

        end
        
        function mask = provide_CONST(mask)
            
        end
        
        function mask = finalize_init(mask)
            temp = load([mask.PARA.point_file_path mask.PARA.filename]);
            mask.PARA.target_lat = temp.target_lat;
            mask.PARA.target_lon = temp.target_lon;
            mask.PARA.target_mask = temp.target_mask;
            
            mask.PARA.target_lat(find(mask.PARA.target_mask(:,1) == 0), :) = []; 
            mask.PARA.target_lon(find(mask.PARA.target_mask(:,1) == 0), :) = [];
            mask.PARA.target_lat_max = nanmax(mask.PARA.target_lat);
            mask.PARA.target_lat_min = nanmin(mask.PARA.target_lat);
            mask.PARA.target_lon_max = nanmax(mask.PARA.target_lon);
            mask.PARA.target_lon_min = nanmin(mask.PARA.target_lon);
        end
        

        function mask = apply_mask(mask)
            
            %[clip_lon, clip_lat] = read_kml(mask, [mask.PARA.kml_file_path mask.PARA.filename]);
            %mask_temp = inpolygon(mask.PARENT.STATVAR.latitude,mask.PARENT.STATVAR.longitude, clip_lat,clip_lon);
            
            
            mask_temp = mask.PARENT.STATVAR.mask .*0;
            %proximity_score = mask.PARENT.STATVAR.latitude.*0 + 1e20;
            if ~(mask.PARA.target_lat_max < nanmin(mask.PARENT.STATVAR.latitude(:))) && ~(mask.PARA.target_lat_min > nanmax(mask.PARENT.STATVAR.latitude(:))) ...
                            && ~(mask.PARA.target_lon_max < nanmin(mask.PARENT.STATVAR.longitude(:))) && ~(mask.PARA.target_lon_min > nanmax(mask.PARENT.STATVAR.longitude(:)))
                        
                for i=1:size(mask.PARA.target_lat,1)
                    if mask.PARA.target_lat(i) > nanmin(mask.PARENT.STATVAR.latitude(:))-0.1 && mask.PARA.target_lat(i) < nanmax(mask.PARENT.STATVAR.latitude(:))+ 0.1 ...
                            && mask.PARA.target_lon(i) > nanmin(mask.PARENT.STATVAR.longitude(:)) - 0.1 && mask.PARA.target_lon(i) < nanmax(mask.PARENT.STATVAR.longitude(:)) + 0.1
                        
                        score=(mask.PARA.target_lon(i) - mask.PARENT.STATVAR.longitude).^2 + (mask.PARA.target_lat(i) - mask.PARENT.STATVAR.latitude).^2;
                        [mini, posi] = min(score(:));
                                                
                        if mini < 0.02
                            mask_temp(posi) = 1;
                            mask.PARENT.STATVAR.properties = [mask.PARENT.STATVAR.properties; [i mini mask.PARENT.PARA.horizontal mask.PARENT.PARA.vertical mask.PARENT.STATVAR.key(posi) ] ];
                        end
                    end
                end
            end
            
            if mask.PARA.additive
                mask.PARENT.STATVAR.mask = mask.PARENT.STATVAR.mask | mask_temp;
            else
                mask.PARENT.STATVAR.mask = mask.PARENT.STATVAR.mask & mask_temp;
            end

        end
        
        
        
        %-------------param file generation-----
        function mask = param_file_info(mask)
            mask = provide_PARA(mask);
            
            mask.PARA.STATVAR = [];
            mask.PARA.class_category = 'MASK';
            mask.PARA.options = [];

            mask.PARA.comment.point_file_path = {'folder where mat-file with target points '};
            mask.PARA.comment.filename = {'name of mat-file with targte points'};
            
            mask.PARA.comment.additive = {'1: region inside kml track added to existing selection; 0: region inside kml track subtracted from existing selection'};
            mask.PARA.default_value.additive = {0};
        end

    end
end

