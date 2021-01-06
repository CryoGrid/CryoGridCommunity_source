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
            run_info = copy(provider.RUN_INFO_CLASS);
            run_info.PPROVIDER = provider;
        end

    end
end

