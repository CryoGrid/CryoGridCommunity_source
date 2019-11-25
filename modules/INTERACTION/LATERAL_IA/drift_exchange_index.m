function drift_index = drift_exchange_index(lateral)
% negative index  - fraction of erodable snow removed
% index =  0 - no snow eroded/deposited
% positive index - fraction of drifting snow which is deposited
    status = lateral.STATUS.snow;
    exposures = lateral.TEMP.exposures;
    area = lateral.PARA.area;
    delta = lateral.PARA.delta;
    drift_index = zeros(1,numlabs);
    
    for i = 1:length(area)
        above = find(exposures > exposures(i) + delta);
        below = find(exposures < exposures(i) - delta);
        if isempty(above) && isempty(below)
            drift_index(i) = 0;
        else
            drift_index(i) = (sum(area(above)) - sum(area(below)))/(sum(area(above)) + sum(area(below)));
        end
    end
    drift_index(drift_index>0) = drift_index(drift_index>0)./sum(drift_index(drift_index>0));
    drift_index = drift_index.*status;
end
