run_list = {'submarine_Benchmark_Testlocation1', 'submarine_Benchmark_Testlocation2', 'submarine_Benchmark_Testlocation3', ...
    'submarine_Benchmark_Testlocation4'};

% run model, looping over run_list
delete(gcp('nocreate')) % kill pre-existing jobs
parpool(min(length(run_list), 8))
parfor i = 1:length(run_list)
    run_number = run_list{i};
    fprintf('Running >>%s<< \n', run_number)
    ok = CG_submarine_main(run_number);
    
    %copy files into zip folder
    try
    resfolder = strcat('results', filesep, run_number, filesep, 'Code.zip');
    zip(resfolder,{'CG_submarine_main.m','CG_submarine_main.m', 'display_ground_SubseaPFwSalt.m', ...
                        strcat('modules', filesep, '*')}); %strcat('forcing', filesep, '*')});
                    
    mkdir(strcat(filesep, 'geo5', filesep, 'HGF_JRG100', filesep, 'CryoGridBenchmarkTests', filesep, 'modularv02', filesep, run_number))
    copyfile(strcat('results', filesep, run_number, filesep, '*'),...
        strcat(filesep, 'geo5', filesep, 'HGF_JRG100', filesep, 'CryoGridBenchmarkTests', filesep, 'modularv02', filesep, run_number))
    catch
        keyboard
        fprintf('could not save files \n')
    end
    
end
