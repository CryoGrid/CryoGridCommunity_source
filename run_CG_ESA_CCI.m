modules_path = 'modules';
addpath(genpath(modules_path));


init_format = 'EXCEL'; %EXCEL or YAML
run_name = 'ESA_CCI'; %parameter file name and result directory 
constant_file = 'CONSTANTS_excel'; %file with constants
result_path = '../results/';  %with trailing backslash
forcing_path = fullfile ('./forcing/');


provider = PROVIDER;
provider = assign_paths(provider, init_format, run_name, result_path, constant_file, forcing_path);
provider = read_const(provider);
provider = read_parameters(provider);


%create the RUN_INFO class
[run_info, provider] = run_model(provider);

run_info = finalize_init(run_info);

%[run_info, tile] = run_preproc(run_info);



%[run_info, tile] = run_model(run_info);

  number_of_tiles = ceil(run_info.PARA.total_number_of_cells ./ run_info.PARA.number_of_cells_per_tile);
            
  for run_index = 1:number_of_tiles
      
      disp(['running range index ' num2str(run_index)])
      
      range = [(run_index-1).*run_info.PARA.number_of_cells_per_tile+1:min(run_index.*run_info.PARA.number_of_cells_per_tile, run_info.PARA.total_number_of_cells)]';
      
      for i=1:size(run_info.PARA.tile_class,1)
          disp(['running tile number ' num2str(i)])
          for j=1:run_info.PARA.number_of_runs(i,1)
              disp(['running round ' num2str(j)])
              
              %load the next tile from the PROVIDER
              tile = copy(run_info.PPROVIDER.CLASSES.(run_info.PARA.tile_class{i,1}){run_info.PARA.tile_class_index(i,1),1});
              tile.PARA.number_of_realizations = size(range,1);
              tile.PARA.range = range;
              
              tile.PARA.geothermal = run_info.STATVAR.geothermal(range,1);
              
              
              %REMOVE
              %                         [~, pos] = max(run_info.STATVAR.landcover(range,1),[], 2);
              %                         tile.PARA.stratigraphy = pos(1,1) .*0 +1;
              %REMOVE
              
              %tile = copy(run_info.PPROVIDER.CLASSES.(run_info.PARA.tile_class){run_info.PARA.tile_class_index,1});
              tile.RUN_INFO = run_info;
              
              tile = finalize_init(tile); %here, tile can still access a potentially existing tile through til.RUN_INFO.TILE
              
              tile = run_model(tile);


          end
      end
  end




%   rest is equivalent to 

%   for RUN_INFO class RUN_1D_STANDARD

% %create the TILE
% tile = copy(run_info.PPROVIDER.CLASSES.(run_info.PARA.tile_class){run_info.PARA.tile_class_index,1});
% tile.RUN_INFO = run_info;
% run_info.TILE = tile;



%create the stratigraphy
%tile = finalize_init(tile);

%    rest is equivalent to 
%tile = run_model(tile);
%   for TILE class TILE_1D_STANDARD




% %=========================================================================
% %TIME INTEGRATION
% %=========================================================================
% while tile.t < tile.FORCING.PARA.end_time
%         
%     %interpolate focing data to time t
%     tile = interpolate_forcing_tile(tile);
% 
%     %upper boundar condition (uppermost class only)
%     TOP.NEXT = get_boundary_condition_u(TOP.NEXT, tile);
%     
%     %set fluxes between classes in the stratigrapht
%     CURRENT = TOP.NEXT;
%     while ~isequal(CURRENT.NEXT, BOTTOM)
%         get_boundary_condition_m(CURRENT.IA_NEXT, tile); %call interaction class function
%         CURRENT = CURRENT.NEXT;
%     end
% 
%     %lower boundary condition (lowermost class)
%     CURRENT = get_boundary_condition_l(CURRENT,  tile);  %At this point, CURRENT is equal to BOTTOM_CLASS
% 
%     %calculate spatial derivatives
%     CURRENT = TOP.NEXT;
%     while ~isequal(CURRENT, BOTTOM)
%         CURRENT = get_derivatives_prognostic(CURRENT, tile);
%         CURRENT = CURRENT.NEXT;
%     end
%     
%     %calculate timestep [second]
%     CURRENT = TOP.NEXT;
%     tile.timestep = 1e8;
%     while ~isequal(CURRENT, BOTTOM)
%         tile.timestep = min(tile.timestep, get_timestep(CURRENT, tile));
%         CURRENT = CURRENT.NEXT;
%     end
%     %tile.timestep = min(4, tile.timestep);
%     tile.next_break_time = min(tile.LATERAL.IA_TIME, tile.OUT.OUTPUT_TIME);
%     tile.timestep = min(tile.timestep, (tile.next_break_time - tile.t).*tile.CONST.day_sec);
%     
%     %prognostic step - integrate prognostic variables in time
%     CURRENT = TOP.NEXT;
%     while ~isequal(CURRENT, BOTTOM)
%         CURRENT = advance_prognostic(CURRENT, tile);
%         CURRENT = CURRENT.NEXT;
%     end
% %     disp(TOP_CLASS.TEMP.d_water(1,1))
% %test=[test; [TOP_CLASS.STATVAR.T(1,1) TOP_CLASS.STATVAR.waterIce(1,1) TOP_CLASS.STATVAR.waterPotential(1,1) TOP_CLASS.TEMP.d_water(1,1)] tile.timestep];
%            
%     %diagnostic step - compute diagnostic variables
%     TOP.NEXT = compute_diagnostic_first_cell(TOP.NEXT, tile); %calculate Lstar, only uppermost class
%     CURRENT = BOTTOM.PREVIOUS;
%     while ~isequal(CURRENT, TOP)
%         CURRENT = compute_diagnostic(CURRENT, tile);
%         CURRENT = CURRENT.PREVIOUS;
%     end
%     
%     %triggers
%     CURRENT = TOP.NEXT;
%     while ~isequal(CURRENT, BOTTOM)
%         CURRENT = check_trigger(CURRENT, tile);
%         CURRENT = CURRENT.NEXT;
%     end
%     
%     tile = interact_lateral(tile);
%     
%     %set TOP_CLASS and BOTTOM_CLASS for convenient access
%     TOP_CLASS = TOP.NEXT; 
%     BOTTOM_CLASS = BOTTOM.PREVIOUS;
%     
%     
% 
%     %update time variable t
%     tile.t = tile.t + tile.timestep./tile.CONST.day_sec;
%     
%     %model 
%     tile = store_OUT_tile(tile);
%    
% end
% 
