
classdef TILE_1D_standard_overwrite_OUT_FORCING < TILE_1D_standard
    
    
    
    methods
        
       

        function tile = provide_PARA(tile)

            tile.PARA.forcing_class = [];
            tile.PARA.forcing_class_index = [];
            tile.PARA.out_class = [];
            tile.PARA.out_class_index = [];

        end
        

       
        
        %assemble the stratigraphy
        function tile = finalize_init(tile)        
            PARA2 = tile.PARA;
            fields = fieldnames(tile.RUN_INFO.TILE);
            for i=1:size(fields,1)
                if ~strcmp(fields{i,1}, 'OUT') || ~strcmp(fields{i,1}, 'FORCING') 
                    tile.(fields{i,1}) = tile.RUN_INFO.TILE.(fields{i,1});
                end
            end
            fields = fieldnames(PARA2);
            for i=1:size(fields,1)
                tile.PARA.(fields{i,1}) = PARA2.(fields{i,1});
            end
            
            %1. forcing
            tile.FORCING = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.forcing_class){tile.PARA.forcing_class_index,1});
            tile.FORCING = finalize_init(tile.FORCING, tile);

            %10. assign time, etc.
            tile.t = tile.FORCING.PARA.start_time;

            %12. assign OUT classes
            tile.OUT = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.out_class){tile.PARA.out_class_index,1});
            tile.OUT = finalize_init(tile.OUT, tile);

        end
        

  
        

    end
end



