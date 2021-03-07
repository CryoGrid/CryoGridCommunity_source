modules_path = 'modules';
addpath(genpath(modules_path));

%-----------------------------
% modified by user
init_format = 'EXCEL3D'; %EXCEL or YAML
%init_format = 'EXCEL';
%run_name = 'Paiku'; %parameter file name and result directory 
%run_name = 'Herschell_test';
%run_name = 'ExperimentHansen2004';
%run_name = 'test_salt';
%run_name = 'revision_paper_juditha';
%run_name = 'example3';
%run_name = 'test_forcing_spinup';
run_name = 'test_BGC_3';
constant_file = 'CONSTANTS_excel'; %file with constants
result_path = '../results/';  %with trailing backslash
%forcing_path = fullfile ('./forcing/');
% end modified by user
%------------------------

%providers
provider = PROVIDER;
provider = assign_paths(provider, init_format, run_name, result_path, constant_file);
provider = read_const(provider);
provider = read_parameters(provider);


% %creates the RUN_INFO class
 [run_info, provider] = run_model(provider);

 [run_info, tile] = run_model(run_info);


