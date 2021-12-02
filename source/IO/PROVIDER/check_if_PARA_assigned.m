function check_if_PARA_assigned(CG_class)

if ~isempty(CG_class.PARA)
    parameters = fieldnames(CG_class.PARA);
    for i=1:size(parameters,1)
        if size(CG_class.PARA.(parameters{i,1}),1) == 0 && size(CG_class.PARA.(parameters{i,1}),2) == 0
            disp(['Warning: PARA ' parameters{i,1} ' in class ' class(CG_class) ' not assigned'])
        end
    end
end