function ia_class = get_IA_class(class1, class2) %change this later, class1 and 2 are strings
 
 
%ia_class = IA_HEAT_SALT;
 
%replace by matrix, this is stupid
 
if strcmp((class1), 'GROUND_freeW_seb') && strcmp((class2), 'GROUND_freeW_seb')
    ia_class = IA_HEAT();
elseif strcmp((class1), 'GROUND_freeW_seb_snow') && strcmp((class2), 'GROUND_freeW_seb_snow')
    ia_class = IA_HEAT();    
elseif (strcmp((class1), 'GROUND_freeW_seb_snow') && strcmp((class2), 'GROUND_freeW_seb')) || (strcmp((class2), 'GROUND_freeW_seb_snow') && strcmp((class1), 'GROUND_freeW_seb'))
    ia_class = IA_HEAT();
elseif (strcmp((class1), 'GROUND_freeW_seb_snow') && strcmp((class2), 'SNOW_seb_simple')) || (strcmp((class2), 'GROUND_freeW_seb_snow') && strcmp((class1), 'SNOW_seb_simple'))
    ia_class = IA_HEAT();
elseif (strcmp((class1), 'GROUND_freeW_seb') && strcmp((class2), 'SNOW_seb_simple')) || (strcmp((class2), 'GROUND_freeW_seb') && strcmp((class1), 'SNOW_seb_simple'))
    ia_class = IA_HEAT();
elseif strcmp((class1), 'GROUND_fcSimple_salt_seb_snow') && strcmp((class2), 'SNOW_seb_simple') || strcmp((class2), 'GROUND_fcSimple_salt_seb_snow') && strcmp((class1), 'SNOW_seb_simple')
    ia_class = IA_HEAT_SALT();
elseif strcmp((class1), 'GROUND_fcSimple_salt_seb_snow') && strcmp((class2), 'GROUND_freeW_seb') || strcmp((class2), 'GROUND_fcSimple_salt_seb_snow') && strcmp((class1), 'GROUND_freeW_seb')
    ia_class = IA_HEAT_SALT();
elseif strcmp((class1), 'GROUND_fcSimple_salt_seb') && strcmp((class2), 'GROUND_freeW_seb') || strcmp((class2), 'GROUND_fcSimple_salt_seb') && strcmp((class1), 'GROUND_freeW_seb')
    ia_class = IA_HEAT_SALT();
elseif strcmp((class1), 'GROUND_fcSimple_salt_seb') && strcmp((class2), 'GROUND_freeW_seb_snow') || strcmp((class2), 'GROUND_fcSimple_salt_seb') && strcmp((class1), 'GROUND_freeW_seb_snow')
    ia_class = IA_HEAT_SALT();
else
    ia_class=0;
    disp('combination of classes not supported')
end
end