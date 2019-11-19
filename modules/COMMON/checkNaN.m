function ground = checkNaN(ground)

names = fieldnames(ground.TEMP);

    for i = 1:size(names,1)
        if isnan(ground.TEMP.(names{i,1}))
            keyboard
        elseif isinf(ground.TEMP.(names{i,1}))
            keyboard
        end
    end

names = fieldnames(ground.STATVAR);

    for i = 1:size(names,1)
        if isnan(ground.STATVAR.(names{i,1}))
            keyboard
        elseif min(ground.STATVAR.T) < -50 
            keyboard
        end
    end

end