function check_if_CONST_assigned(CG_class)

if ~isempty(CG_class.CONST)
    parameters = fieldnames(CG_class.CONST);
    for i=1:size(parameters,1)
        if size(CG_class.CONST.(parameters{i,1}),1) == 0 && size(CG_class.CONST.(parameters{i,1}),2) == 0
            disp(['Warning: CONST ' parameters{i,1} ' in class ' class(CG_class) ' not assigned'])
        end
    end
    
end