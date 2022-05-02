%========================================================================
% CryoGrid LATERAL_IA class LAT_SEEPAGE_FACE_WATER 
% LAT_SEEPAGE_FACE_WATER simulates lateral water flow through a seepage 
% face with defined upper and lower elevations (absolute elevation, 
% not relative to the surface!), as well as contact lengths (i.e. width) 
% and distance from the subsurface class.
% At the seepage face, air pressure is assumed, while the water pressure 
% in the stratigraphy corepsonds to the hydrostatic pressure, i.e. water
% flow only occurs from a saturated zone.
% S. Westermann, Oct 2020
%========================================================================


classdef LAT_SEEPAGE_FACE_WATER < BASE_LATERAL

    
    methods
        
        %----mandatory functions---------------
        %----initialization--------------------

        function lateral = provide_CONST(lateral)
            lateral.CONST.day_sec = [];
        end
        
        function lateral = provide_PARA(lateral)
            lateral.PARA.upperElevation = []; %Inf;
            lateral.PARA.lowerElevation = []; %20;
            lateral.PARA.hardBottom_cutoff = []; %0.03; %hard bottom if saturated and water content below
            lateral.PARA.distance_seepageFace = []; %1; 
            lateral.PARA.seepage_contact_length = []; %4;
        end
        
        function lateral = provide_STATVAR(lateral)
            lateral.STATVAR.subsurface_run_off = [];
        end

        function lateral = finalize_init(lateral, tile)
            lateral.STATVAR.subsurface_run_off = 0;
        end

        %----time integration------
        
        %only push function needed
        function lateral = push(lateral, tile)
            
            CURRENT = lateral.PARENT.TOP.NEXT;
            lateral.TEMP.head = 0;
            while ~(strcmp(class(CURRENT), 'Bottom'))
                CURRENT = lateral_push_remove_water_seepage(CURRENT, lateral);
                CURRENT = compute_diagnostic(CURRENT, tile);
                CURRENT = CURRENT.NEXT;
            end
            
        end
        
        
        function lateral = set_ACTIVE(lateral, i, t)
            lateral.PARENT.ACTIVE(i,1) = 1;
        end
        
        function lateral = get_derivatives(lateral, tile)
            
        end
        
        function lateral = pull(lateral, tile)
            
        end
        
        
        %-------------param file generation-----
        function ground = param_file_info(ground)
            ground = param_file_info@BASE_LATERAL(ground);
            
            ground.PARA.class_category = 'LATERAL_IA';
            
            ground.PARA.options = [];
            ground.PARA.STATVAR = [];
            
            ground.PARA.default_value.upperElevation = {10000};
            ground.PARA.comment.upperElevation = {'upper elevation of seepage face [meter a.s.l.]'};
            
            ground.PARA.default_value.lowerElevation = {0};
            ground.PARA.comment.lowerElevation = {'lower elevation of seepage face [meter a.s.l.]'};
            
            ground.PARA.default_value.hardBottom_cutoff = {0.03};
            ground.PARA.comment.hardBottom_cutoff = {'hard bottom  = no water flow if saturated and water content below [vol. water content, -]'};
            
            ground.PARA.default_value.distance_seepageFace = {1};
            ground.PARA.comment.distance_seepageFace ={'distance to seepage face [m]'};
            
            ground.PARA.default_value.seepage_contact_length = {1};
            ground.PARA.comment.seepage_contact_length = {'lateral contact length (=width) of seepage face [m]'};
        end
    end
    
end


