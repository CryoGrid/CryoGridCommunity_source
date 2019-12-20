function lateral = precondition_lateral_exchange(lateral,top_class)
CLASS = class(top_class);
if strcmp(CLASS(1:6),'GROUND') && top_class.STATVAR.T(1) <= 0
    precond = 0;
else
    precond = 1;
end
lateral.STATUS.water(labindex) = precond;

if lateral.PARA.ghost == 1
    lateral.STATUS.water(end) = 1;
end

for j = 1:numlabs
    if j ~= labindex
        labSend(class(top_class),j,102);
        labSend(precond,j,103);
    end
end
for j = 1:numlabs
    if j ~= labindex
        class_j = labReceive(j,102);
        precond_j = labReceive(j,103);
        precond_j = precond_j*precond; 
        if strcmp(class_j,CLASS)
            lateral.STATUS.water(j) = precond_j;
        else
            lateral.STATUS.water(j) = 0;
        end
    end
end

end
