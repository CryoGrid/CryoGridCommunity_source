modules_path = 'modules';
addpath(genpath(modules_path));

%-----------------------------
% modified by user
init_format = 'EXCEL'; %EXCEL or YAML
run_name = 'test'; %parameter file name and result directory 
%run_name = 'Herschell_test';
constant_file = 'CONSTANTS_excel'; %file with constants
result_path = '../results/';  %with trailing backslash
forcing_path = fullfile ('./forcing/');
% end modified by user
%------------------------

%providers
provider = PROVIDER;
provider = assign_paths(provider, init_format, run_name, result_path, constant_file, forcing_path);
provider = read_const(provider);
provider = read_parameters(provider);


% %creates the RUN_INFO class
% [run_info, provider] = run_model(provider);
% 
% [run_info, tile] = run_model(run_info);


