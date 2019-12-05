function ground = add_CHILD_snow(ground, class_list, stratigraphy_list)

for i=1:size(stratigraphy_list,1)  %find STRAT_class in the list
    if stratigraphy_list{i,2}==1
        if strcmp(class(stratigraphy_list{i,1}), 'STRAT_classes')
            class_stratigraphy = stratigraphy_list{i,1};
            if isempty(class_stratigraphy.snow_class)
                return
            end
        end
    end
end

for i=1:size(class_list,1)  %find snow in the class list
    if strcmp(class(class_list{i,1}), class_stratigraphy.snow_class.name) && class_list{i,2} == class_stratigraphy.snow_class.index
        snow=copy(class_list{i,1});
    end
end

%replace by matrix

if strcmp(class(ground), 'GROUND_freeW_seb_snow') && strcmp(class(snow), 'SNOW_simple_seb_bucketW')
    ground.IA_CHILD = IA_SNOW_GROUND();
elseif strcmp(class(ground), 'GROUND_freeW_seb_snow') && strcmp(class(snow), 'SNOW_simple_seb_crocus')
    ground.IA_CHILD = IA_SNOW_GROUND_crocus();
elseif strcmp(class(ground), 'GROUND_fcSimple_salt_seb_snow') && strcmp(class(snow), 'SNOW_simple_seb_bucketW')
    ground.IA_CHILD = IA_SNOW_GROUND_fcSimple_salt();
elseif strcmp(class(ground), 'GROUND_fcSimple_salt_seb_snow') && strcmp(class(snow), 'SNOW_crocus_no_inheritance')
    ground.IA_CHILD = IA_SNOW_GROUND_fcSimple_salt_crocus();
elseif strcmp(class(ground), 'GROUND_freeW_bucketW_seb_snow') && strcmp(class(snow), 'SNOW_simple_seb_bucketW')
    ground.IA_CHILD = IA_SNOW_GROUND();
elseif strcmp(class(ground), 'GROUND_freeW_bucketW_seb_snow') && strcmp(class(snow), 'SNOW_simple_seb_crocus')
    ground.IA_CHILD = IA_SNOW_GROUND_crocus();
elseif strcmp(class(ground), 'GROUND_freeW_bucketW_seb_snow') && strcmp(class(snow), 'SNOW_crocus_no_inheritance')
    ground.IA_CHILD = IA_SNOW_GROUND_crocus();
elseif strcmp(class(ground), 'GROUND_freeW_seb_snow') && strcmp(class(snow), 'SNOW_simple_seb')
    ground.IA_CHILD = IA_SNOW_GROUND();
elseif strcmp(class(ground), 'GROUND_freezeC_bucketW_seb_snow') && strcmp(class(snow), 'SNOW_simple_seb_crocus')
    ground.IA_CHILD = IA_SNOW_GROUND_crocus_SW();
elseif strcmp(class(ground), 'GROUND_freezeC_bucketW_seb_snow') && strcmp(class(snow), 'SNOW_simple_seb')
    ground.IA_CHILD = IA_SNOW_GROUND_crocus_SW();
elseif strcmp(class(ground), 'GROUND_fcSimple_salt_seb_snow_test') && strcmp(class(snow), 'SNOW_simple_seb_crocus')
    ground.IA_CHILD = IA_SNOW_GROUND_crocus_SW();
elseif strcmp(class(ground), 'GROUND_fcSimple_salt_seb_snow_test') && strcmp(class(snow), 'SNOW_simple_seb')
    ground.IA_CHILD = IA_SNOW_GROUND_crocus_SW();    

end


CURRENT = ground.IA_CHILD; %change to interaction class
CURRENT.STATUS = 0; %snow initially inactive
CURRENT.IA_PARENT_GROUND = ground;
CURRENT.IA_CHILD_SNOW = snow;
CURRENT.IA_CHILD_SNOW.IA_PARENT = CURRENT;
CURRENT.IA_CHILD_SNOW = initialize_zero_snow(CURRENT.IA_CHILD_SNOW, CURRENT.IA_PARENT_GROUND);

ground.IA_CHILD = CURRENT; %reassign to ground







