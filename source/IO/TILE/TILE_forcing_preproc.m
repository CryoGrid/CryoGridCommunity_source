% base class build a model tile

classdef TILE_forcing_preproc < matlab.mixin.Copyable
    
    properties
        

        BUILDER
        PARA
        RUN_INFO
        FORCING
        CONST
        GRID
        OUT        

    end

    methods

        function tile = provide_PARA(tile)
            
            tile.PARA.forcing_path = [];
            tile.PARA.forcing_file_name = [];
            tile.PARA.start_date = [];
            tile.PARA.end_date = [];
            tile.PARA.number_of_repetitions = [];
            tile.PARA.real_time = []; %repetetion number, 0 otherwise

            tile.PARA.save_path = [];
            tile.PARA.save_file_name = [];
            
        end
        
        function tile = provide_CONST(tile)

        end
        
        function tile = provide_STATVAR(tile)
            
            
        end
        

        %assemble the stratigraphy
        function tile = finalize_init(tile)
            start_date = [];
            end_date = [];
            for i=1:size(tile.PARA.start_date,1)

                if ~isempty(tile.PARA.start_date{i,1}) && ~sum(double(isnan(tile.PARA.start_date{i,1}))>0)
                    start_date = [start_date; datenum(tile.PARA.start_date{i,1}, 'dd.mm.yyyy')];
                else
                    start_date = [start_date; NaN];
                end
                if ~isempty(tile.PARA.end_date{i,1}) && ~sum(double(isnan(tile.PARA.end_date{i,1}))>0)
                    end_date = [end_date; datenum(tile.PARA.end_date{i,1}, 'dd.mm.yyyy')];
                else    
                    end_date = [end_date; NaN];
                end
            end
   
            tile.PARA.start_date = start_date;
            tile.PARA.end_date = end_date;
        end

        function tile = interpolate_forcing_tile(tile)

        end

        function tile = interact_lateral(tile)

        end
        
        function tile = store_OUT_tile(tile)

        end        
        
        
        
        function tile = run_model(tile)
             real_time_index = find(tile.PARA.real_time(:,1) >0 );
             FORCING = load([tile.PARA.forcing_path tile.PARA.forcing_file_name{real_time_index,1}]);
             FORCING = FORCING.FORCING;

             variables = fieldnames(FORCING.data);
             
             if isnan(tile.PARA.start_date(real_time_index,1))
                 start_pos = 1;
             else
                 [~, start_pos] = min(abs(FORCING.data.t_span-tile.PARA.start_date(real_time_index,1)));
             end

             if isnan(tile.PARA.end_date(real_time_index,1))
                 end_pos = size(FORCING.data.t_span,1);
             else
                 [~, end_pos] = min(abs(FORCING.data.t_span-tile.PARA.end_date(real_time_index,1)));
                 end_pos = end_pos-1;
             end
             
             for i=1:size(variables,1)

                 if size(FORCING.data.(variables{i,1}),1)>=end_pos
                     tile.FORCING.data.(variables{i,1}) = FORCING.data.(variables{i,1})(start_pos:end_pos,1);
                 end
             end
            
             %same file backwards
             variables = fieldnames(tile.FORCING.data);
             for ii = tile.PARA.real_time(real_time_index,1)-1:-1:1
                              
                 for i=1:size(variables,1)
                     if ~strcmp(variables{i,1} ,'t_span')
                         tile.FORCING.data.(variables{i,1}) = [FORCING.data.(variables{i,1})(start_pos:end_pos,1); tile.FORCING.data.(variables{i,1})];
                     else
                         offset = FORCING.data.(variables{i,1})(end_pos,1) - tile.FORCING.data.(variables{i,1})(1,1) +(tile.FORCING.data.(variables{i,1})(2,1)-tile.FORCING.data.(variables{i,1})(1,1));
                         tile.FORCING.data.(variables{i,1}) = [FORCING.data.(variables{i,1})(start_pos:end_pos,1)-offset; tile.FORCING.data.(variables{i,1})];
                     end
                 end
             end
             %same file forward
             for ii = tile.PARA.real_time(real_time_index,1)+1:tile.PARA.number_of_repetitions(real_time_index,1)
                 
                 for i=1:size(variables,1)
                     if ~strcmp(variables{i,1} ,'t_span')
                         tile.FORCING.data.(variables{i,1}) = [ tile.FORCING.data.(variables{i,1}); FORCING.data.(variables{i,1})(start_pos:end_pos,1)];
                     else
                         offset = FORCING.data.(variables{i,1})(start_pos,1) - tile.FORCING.data.(variables{i,1})(end,1) - (tile.FORCING.data.(variables{i,1})(2,1)-tile.FORCING.data.(variables{i,1})(1,1));
                         tile.FORCING.data.(variables{i,1}) = [ tile.FORCING.data.(variables{i,1}); FORCING.data.(variables{i,1})(start_pos:end_pos,1)-offset];
                     end
                 end
             end
             
             
             %previous files backwards
             for jj=real_time_index-1:-1:1
                 for ii = 1:tile.PARA.number_of_repetitions(jj,1)
                     FORCING = load([tile.PARA.forcing_path tile.PARA.forcing_file_name{jj,1}]);
                     FORCING = FORCING.FORCING;
                     
                     if isnan(tile.PARA.start_date(jj,1))
                         start_pos = 1;
                     else
                         [~, start_pos] = min(abs(FORCING.data.t_span-tile.PARA.start_date(jj,1)));
                     end
                     
                     if isnan(tile.PARA.end_date(jj,1))
                         end_pos = size(FORCING.data.t_span,1);
                     else
                         [~, end_pos] = min(abs(FORCING.data.t_span-tile.PARA.end_date(jj,1)));
                         end_pos = end_pos-1;
                     end
                     
                     for i=1:size(variables,1)
                     if ~strcmp(variables{i,1} ,'t_span')
                         tile.FORCING.data.(variables{i,1}) = [FORCING.data.(variables{i,1})(start_pos:end_pos,1); tile.FORCING.data.(variables{i,1})];
                     else
                         offset = FORCING.data.(variables{i,1})(end_pos,1) - tile.FORCING.data.(variables{i,1})(1,1) +(tile.FORCING.data.(variables{i,1})(2,1)-tile.FORCING.data.(variables{i,1})(1,1));
                         tile.FORCING.data.(variables{i,1}) = [FORCING.data.(variables{i,1})(start_pos:end_pos,1)-offset; tile.FORCING.data.(variables{i,1})];
                     end
                 end

                 end
             end
             
             
             %next files forward
             for jj=real_time_index+1:size(tile.PARA.real_time,1)
                 for ii = 1:tile.PARA.number_of_repetitions(jj,1)
                     
                     FORCING = load([tile.PARA.forcing_path tile.PARA.forcing_file_name{jj,1}]);
                     FORCING = FORCING.FORCING;
                     
                     if isnan(tile.PARA.start_date(jj,1))
                         start_pos = 1;
                     else
                         [~, start_pos] = min(abs(FORCING.data.t_span-tile.PARA.start_date(jj,1)));
                     end
                     
                     if isnan(tile.PARA.end_date(jj,1))
                         end_pos = size(FORCING.data.t_span,1);
                     else
                         [~, end_pos] = min(abs(FORCING.data.t_span-tile.PARA.end_date(jj,1)));
                         end_pos = end_pos-1;
                     end
                     
                     for i=1:size(variables,1)
                         if ~strcmp(variables{i,1} ,'t_span')
                             tile.FORCING.data.(variables{i,1}) = [ tile.FORCING.data.(variables{i,1}); FORCING.data.(variables{i,1})(start_pos:end_pos,1)];
                         else
                             offset = FORCING.data.(variables{i,1})(start_pos,1) - tile.FORCING.data.(variables{i,1})(end,1) - (tile.FORCING.data.(variables{i,1})(2,1)-tile.FORCING.data.(variables{i,1})(1,1));
                             tile.FORCING.data.(variables{i,1}) = [ tile.FORCING.data.(variables{i,1}); FORCING.data.(variables{i,1})(start_pos:end_pos,1)-offset];
                         end
                     end
                 end
             end
            
             
             FORCING = tile.FORCING;
             save([tile.PARA.save_path tile.PARA.save_file_name], 'FORCING')

        end
        
        

    end
end



