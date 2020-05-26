function ia_class = get_IA_class(below_class, above_class)

classes = [ {'GROUND_fcSimple_salt'};
    {'GROUND_fcSimple_salt_seb'};
    {'GROUND_fcSimple_salt_seb_snow'};
    {'GROUND_fcSimple_salt_seb_snow_test'};
    {'GROUND_freeW_bucketW'};
    {'GROUND_freeW_bucketW_seb'};
    {'GROUND_freeW_bucketW_seb_snow'};
    {'GROUND_freeW_seb'};
    {'GROUND_freeW_seb_snow'};
    {'GROUND_freezeC_bucketW_seb_snow'}; % 10
    {'GROUND_freezeC_bucketW_seb'};    
    {'GROUND_freezeC_bucketW_Xice_seb_snow'};
    {'SNOW_simple_seb'};
    {'SNOW_simple_seb_bucketW'};
    {'SNOW_simple_seb_crocus'}; % 15
    {'SNOW_crocus_no_inheritance'}
    {'GROUND_vegetation'}; % 17    
    {'GROUND_vegetation_snow'};  
    {'SNOW_simple_seb_vegetation'};
    ];

ia=[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;0 0 0 0 0 0 0 2 2 0 0 0 0 0 0 0 0 0 0;0 0 0 0 0 0 0 2 0 0 0 0 0 1 0 1 0 0 0;0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 0 0 0;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;0 0 0 0 0 0 0 3 0 0 0 0 0 0 0 0 0 0 0;0 0 0 0 0 0 0 3 0 0 0 0 0 0 0 0 0 0 0;0 2 2 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0;0 2 0 0 0 0 0 1 1 0 0 0 1 1 1 0 0 0 0;0 0 0 0 0 0 0 0 0 0 0 0 1 1 7 7 6 0 1;0 0 0 0 0 0 0 0 0 0 0 0 0 0 7 0 6 6 0;0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 0 0 0;0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0;0 0 1 0 0 0 4 0 1 0 0 0 0 0 0 0 0 0 0;0 0 0 0 0 0 4 0 1 7 1 0 0 0 0 0 0 1 0;0 0 1 0 0 0 4 0 0 0 0 0 0 0 0 0 0 0 0;0 0 0 0 0 0 0 0 0 6 6 0 0 0 0 0 6 0 0;0 0 0 0 0 0 0 0 0 0 6 0 0 0 1 0 0 0 0;0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0];
% ia=[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;0 0 0 0 0 0 0 2 2 0 0 0 0 0 0 0 0 0 0;0 0 0 0 0 0 0 2 0 0 0 0 0 1 0 1 0 0 0;0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 0 0 0;0 0 0 0 0 0 0 0 0 4 4 0 0 0 0 0 5 0 0;0 0 0 0 0 0 0 4 0 0 0 0 0 0 0 0 5 0 0;0 0 0 0 0 0 0 4 0 0 0 0 0 0 0 0 5 0 0;0 2 2 0 0 0 0 1 1 0 0 0 0 0 0 0 5 0 0;0 2 0 0 4 0 0 4 4 0 0 0 1 1 1 0 0 0 0;0 0 0 0 4 0 0 0 0 4 4 0 1 1 1 1 0 1 0;0 0 0 0 4 0 0 0 0 4 4 0 0 0 0 0 5 0 5;0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 0 0 0;0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0;0 0 1 0 0 0 4 0 1 0 0 0 0 0 0 0 0 0 0;0 0 0 0 0 0 4 0 1 0 0 0 0 0 0 0 0 0 0;0 0 1 0 0 0 4 0 0 0 0 0 0 0 0 0 0 0 0;0 0 0 0 5 5 0 5 0 1 5 0 0 0 0 0 5 0 0;0 0 0 0 5 5 0 5 0 1 5 0 0 0 0 0 5 0 1];
 
for i = 1:length(classes)
    if strcmp(above_class,classes{i})
        id_a = i;
    end
    if strcmp(below_class,classes{i})
        id_b = i;
    end
end

if ia(id_b,id_a) == 1
    ia_class = IA_HEAT();
elseif ia(id_b,id_a) == 2
    ia_class = IA_HEAT_SALT();
elseif ia(id_b,id_a) == 3
    ia_class = IA_HEAT_WATER_1();
elseif ia(id_b,id_a) == 4
    ia_class = IA_HEAT_WATER_SNOW_GROUND();
elseif ia(id_b,id_a) == 5
    ia_class = IA_HEAT_VEGETATION();
elseif ia(id_b,id_a) == 6 %CHANGED 20200121
    ia_class = IA_HEAT_GROUND_VEGETATION();
elseif ia(id_b,id_a) == 7 %CHANGED 20200121
    ia_class = IA_HEAT_GROUND_SNOW_VEGETATION();
else
    ia_class=0;
    disp('combination of classes not supported')
end

%{
for j = 1:12
for i = 1:12
ia_class = get_IA_class_old(classes{i},classes{j});
if ia_class == 0
ia(i,j) = 0;
elseif strcmp(class(ia_class),'IA_HEAT')
ia(i,j) = 1;
elseif strcmp(class(ia_class),'IA_HEAT_SALT')
ia(i,j) = 2;
elseif strcmp(class(ia_class),'IA_HEAT_WATER_1')
ia(i,j) = 3;
elseif strcmp(class(ia_class),'IA_HEAT_WATER_SNOW_GROUND')
    ia(i,j) = 4;
end
end
end
%}

end