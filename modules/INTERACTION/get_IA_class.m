function ia_class = get_IA_class(class1, class2) %change this later, class1 and 2 are strings


%ia_class = IA_HEAT_SALT;

%replace by matrix, this is stupid

if strcmp((class1), 'GROUND_freeW_seb') && strcmp((class2), 'GROUND_freeW_seb')
    ia_class = IA_HEAT();
elseif strcmp((class1), 'GROUND_freeW_seb_snow') && strcmp((class2), 'GROUND_freeW_seb_snow')
    ia_class = IA_HEAT();
elseif (strcmp((class1), 'GROUND_freeW_seb_snow') && strcmp((class2), 'GROUND_freeW_seb')) || (strcmp((class2), 'GROUND_freeW_seb_snow') && strcmp((class1), 'GROUND_freeW_seb'))
    ia_class = IA_HEAT();
elseif (strcmp((class1), 'GROUND_freeW_seb_snow') && strcmp((class2), 'SNOW_simple_seb_bucketW')) || (strcmp((class2), 'GROUND_freeW_seb_snow') && strcmp((class1), 'SNOW_simple_seb_bucketW'))
    ia_class = IA_HEAT();
elseif (strcmp((class1), 'GROUND_freeW_seb_snow') && strcmp((class2), 'SNOW_sim66ple_seb_crocus')) || (strcmp((class2), 'GROUND_freeW_seb_snow') && strcmp((class1), 'SNOW_simple_seb_crocus'))
    ia_class = IA_HEAT();
elseif strcmp((class1), 'GROUND_fcSimple_salt_seb_snow') && strcmp((class2), 'SNOW_simple_seb_bucketW') || strcmp((class2), 'GROUND_fcSimple_salt_seb_snow') && strcmp((class1), 'SNOW_simple_seb_bucketW')
    ia_class = IA_HEAT();
elseif strcmp((class1), 'GROUND_fcSimple_salt_seb_snow') && strcmp((class2), 'GROUND_freeW_seb') || strcmp((class2), 'GROUND_fcSimple_salt_seb_snow') && strcmp((class1), 'GROUND_freeW_seb')
    ia_class = IA_HEAT_SALT();
elseif strcmp((class1), 'GROUND_fcSimple_salt_seb') && strcmp((class2), 'GROUND_freeW_seb') || strcmp((class2), 'GROUND_fcSimple_salt_seb') && strcmp((class1), 'GROUND_freeW_seb')
    ia_class = IA_HEAT_SALT();
elseif strcmp((class1), 'GROUND_fcSimple_salt_seb') && strcmp((class2), 'GROUND_freeW_seb_snow') || strcmp((class2), 'GROUND_fcSimple_salt_seb') && strcmp((class1), 'GROUND_freeW_seb_snow')
    ia_class = IA_HEAT_SALT();
elseif strcmp((class1), 'GROUND_freeW_bucketW_seb') && strcmp((class2), 'GROUND_freeW_seb')
    ia_class = IA_HEAT_WATER_1();
elseif strcmp((class1), 'GROUND_freeW_bucketW_seb_snow') && strcmp((class2), 'GROUND_freeW_seb')
    ia_class = IA_HEAT_WATER_1();
elseif strcmp((class1), 'SNOW_simple_seb_bucketW') && strcmp((class2), 'GROUND_freeW_bucketW_seb_snow')
    ia_class = IA_HEAT_WATER_SNOW_GROUND();
elseif strcmp((class1), 'SNOW_simple_seb_crocus') && strcmp((class2), 'GROUND_freeW_bucketW_seb_snow')
    ia_class = IA_HEAT_WATER_SNOW_GROUND();
elseif (strcmp((class1), 'GROUND_freeW_seb_snow') && strcmp((class2), 'SNOW_simple_seb')) || (strcmp((class2), 'GROUND_freeW_seb_snow') && strcmp((class1), 'SNOW_simple_seb'))
    ia_class = IA_HEAT();
else
    ia_class=0;
    disp('combination of classes not supported')
end
end