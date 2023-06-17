classdef BASE_PROVIDER
    
    properties
        PARA
        CONST
        CLASSES %struct of all classes
        STORAGE
        RUN_INFO_CLASS %RUN_INFO class to start
    end
    
    methods

        function [run_info, provider] = run_model(provider)
            provider = check_PARA_CONST(provider);
            run_info = copy(provider.RUN_INFO_CLASS);
            run_info.PPROVIDER = provider;
            run_info = finalize_init(run_info);
        end

        %check if all PARA and CONST are assigned; generic function used,
        %which can be overwritten in individual classes, like
        %TILE_1D_standard, if PARA itself depends on input PARA
        function provider = check_PARA_CONST(provider)
            class_names = fieldnames(provider.CLASSES);
            for i=1:size(class_names,1)
                for j=1:size(provider.CLASSES.(class_names{i,1}),1)
                    check_if_PARA_assigned(provider.CLASSES.(class_names{i,1}){j,1});
                    check_if_CONST_assigned(provider.CLASSES.(class_names{i,1}){j,1});
                end
            end
        end

    end
end

