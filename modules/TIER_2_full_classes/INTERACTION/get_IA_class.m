function ia_class = get_IA_class(above_class, below_class)

%change Tier 4 cases to their TIER3 and TIER2 base classes 

%[above_class, below_class] = get_IA_class_TIER4(above_class, below_class);


classes = [ {'GROUND_freeW_seb'};
    {'GROUND_freeW_seb_snow'};
    {'GROUND_freeW_bucketW_seb'};
    {'GROUND_freeW_bucketW_seb_snow'};
    {'GROUND_freeW_bucketW_convection_seb'};
    {'GROUND_freeW_bucketW_convection_seb_snow'};
    {'GROUND_freezeC_seb'};
    {'GROUND_freezeC_seb_snow'};
    {'GROUND_freezeC_bucketW_seb'};
    {'GROUND_freezeC_bucketW_seb_snow'};
    {'GROUND_freezeC_bucketW_Xice_seb'};
    {'GROUND_freezeC_bucketW_Xice_seb_snow'};
    {'GROUND_freezeC_RichardsEqW_seb'};
    {'GROUND_freezeC_RichardsEqW_seb_snow'};
    {'GROUND_freezeC_RichardsEqW_Xice_seb'};
    {'GROUND_freezeC_RichardsEqW_Xice_seb_snow'};
    {'GROUND_fcSimple_salt_seb'};
    {'GROUND_fcSimple_salt_seb_snow'};
    {'LAKE_simple_seb'};
    {'LAKE_simple_seb_snow'};
    {'LAKE_simple_unfrozen_seb'};
    {'LAKE_simple_bucketW_seb'};
    {'LAKE_simple_bucketW_seb_snow'};
    {'LAKE_simple_unfrozen_bucketW_seb'};
    {'SNOW_simple_seb'};
    {'SNOW_simple_bucketW_seb'};
    {'SNOW_crocus_bucketW_seb'};
    {'SNOW_crocus2_bucketW_seb'}];



for i = 1:length(classes)
    if strcmp(below_class,classes{i})
        id_a = i;
    end
    if strcmp(above_class,classes{i})
        id_b = i;
    end
end



% ia=[1	0	3	3	3	3	1	0	3	3	3	3	6	6	0	0	0	0	0	0	0	0	0	0;
%     1	1	3	3	3	3	1	1	3	3	3	3	6	6	0	0	0	0	0	0	0	0	0	0;
%     2	0	4	4	4	4	2	0	4	4	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
%     2	0	4	0	4	0	2	0	4	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
%     2	0	4	4	0	0	2	0	4	4	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
%     2	0	4	0	0	0	2	0	4	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
%     1	0	3	3	3	3	1	0	3	3	3	3	6	6	0	0	0	0	0	0	0	0	0	0;
%     1	1	3	3	3	3	1	1	3	3	3	3	6	6	0	0	0	0	0	0	0	0	0	0;
%     2	0	4	4	4	4	2	0	4	4	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
%     2	0	4	0	4	0	2	0	4	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
%     2	0	0	0	0	0	2	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
%     2	0	0	0	0	0	2	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
%     5	0	0	0	0	0	5	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
%     5	0	0	0	0	0	5	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
%     1	1	0	0	0	0	1	1	0	0	0	0	0	0	0	0	101	0	0	0	0	0	0	0;
%     1	1	0	0	0	0	1	1	0	0	0	0	0	0	0	0	101	0	0	0	0	0	0	0;
%     7	7	0	0	0	0	7	7	0	0	0	0	0	0	101	101	0	0	0	0	0	0	0	0;
%     8	8	9	9	9	9	8	8	9	9	13	13	0	0	0	0	0	101	101	101	0	0	0	0;
%     8	8	9	9	9	9	8	8	9	9	13	13	0	0	0	0	0	101	101	101	0	0	0	0;
%     8	8	9	9	9	9	8	8	9	9	13	13	0	0	0	0	0	101	101	101	0	0	0	0;
%     0	1	0	3	0	3	0	1	0	3	0	3	0	1	0	1	0	0	3	0	0	0	0	0;
%     0	2	0	10	0	10	0	2	0	10	0	11	0	2	0	2	0	0	12	0	0	0	0	0;
%     0	2	0	10	0	10	0	2	0	10	0	11	0	2	0	2	0	0	12	0	0	0	0	0;
%     0	0	0	0	0	0	0	0	0	10	0	11	0	0	0	0	0	0	12	0	0	0	0	0];

% ia=[1	0	3	3	3	3	1	0	3	3	3	3	0	0	6	6	0	0	0	0	0	0	0	0	0	0;
%     1	1	3	3	3	3	1	1	3	3	3	3	0	0	6	6	0	0	0	0	0	0	0	0	0	0;
%     2	0	4	4	4	4	2	0	4	4	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
%     2	0	4	0	4	0	2	0	4	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
%     2	0	4	4	0	0	2	0	4	4	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
%     2	0	4	0	0	0	2	0	4	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
%     1	0	3	3	3	3	1	0	3	3	3	3	0	0	6	6	0	0	0	0	0	0	0	0	0	0;
%     1	1	3	3	3	3	1	1	3	3	3	3	0	0	6	6	0	0	0	0	0	0	0	0	0	0;
%     2	0	4	4	4	4	2	0	4	4	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
%     2	0	4	0	4	0	2	0	4	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
%     2	0	0	0	0	0	2	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
%     2	0	0	0	0	0	2	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
%     2	2	0	0	0	0	2	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
%     2	2	0	0	0	0	2	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
%     5	0	0	0	0	0	5	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
%     5	0	0	0	0	0	5	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
%     1	1	0	0	0	0	1	1	0	0	0	0	0	0	0	0	0	0	101	0	0	0	0	0	0	0;
%     1	1	0	0	0	0	1	1	0	0	0	0	0	0	0	0	0	0	101	0	0	0	0	0	0	0;
%     7	7	0	0	0	0	7	7	0	0	0	0	0	0	0	0	101	101	0	0	0	0	0	0	0	0;
%     8	8	9	9	9	9	8	8	9	9	13	13	16	16	0	0	0	0	0	101	101	101	0	0	0	0;
%     8	8	9	9	9	9	8	8	9	9	13	13	16	16	0	0	0	0	0	101	101	101	0	0	0	0;
%     8	8	9	9	9	9	8	8	9	9	13	13	16	16	0	0	0	0	0	101	101	101	0	0	0	0;
%     0	1	0	3	0	3	0	1	0	3	0	3	0	3	0	1	0	1	0	0	3	0	0	0	0	0;
%     0	2	0	10	0	10	0	2	0	10	0	11	0	15	0	2	0	2	0	0	12	0	0	0	0	0;
%     0	2	0	10	0	10	0	2	0	10	0	11	0	15	0	2	0	2	0	0	12	0	0	0	0	0;
%     0	0	0	0	0	0	0	0	0	10	0	11	0	15	0	0	0	0	0	0	12	0	0	0	0	0];

ia=[1	0	3	3	3	3	1	0	3	3	3	3	0	0	0	0	6	6	0	0	0	0	0	0	0	0	0	0;
    1	1	3	3	3	3	1	1	3	3	3	3	0	0	0	0	6	6	0	0	0	0	0	0	0	0	0	0;
    2	0	4	4	4	4	2	0	4	4	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
    2	0	4	0	4	0	2	0	4	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
    2	0	4	4	0	0	2	0	4	4	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
    2	0	4	0	0	0	2	0	4	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
    1	0	3	3	3	3	1	0	3	3	3	3	0	0	0	0	6	6	0	0	0	0	0	0	0	0	0	0;
    1	1	3	3	3	3	1	1	3	3	3	3	0	0	0	0	6	6	0	0	0	0	0	0	0	0	0	0;
    2	0	4	4	4	4	2	0	4	4	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
    2	0	4	0	4	0	2	0	4	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
    2	0	0	0	0	0	2	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
    2	0	0	0	0	0	2	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
    2	2	0	0	0	0	2	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
    2	2	0	0	0	0	2	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
    2	2	0	0	0	0	2	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
    2	2	0	0	0	0	2	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
    5	0	0	0	0	0	5	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
    5	0	0	0	0	0	5	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
    1	1	0	0	0	0	1	1	0	0	0	0	0	0	0	0	0	0	0	0	101	0	0	0	0	0	0	0;
    1	1	0	0	0	0	1	1	0	0	0	0	0	0	0	0	0	0	0	0	101	0	0	0	0	0	0	0;
    7	7	0	0	0	0	7	7	0	0	0	0	0	0	0	0	0	0	101	101	0	0	0	0	0	0	0	0;
    8	8	9	9	9	9	8	8	9	9	13	13	16	16	19	19	0	0	0	0	0	101	101	101	0	0	0	0;
    8	8	9	9	9	9	8	8	9	9	13	13	16	16	19	19	0	0	0	0	0	101	101	101	0	0	0	0;
    8	8	9	9	9	9	8	8	9	9	13	13	16	16	19	19	0	0	0	0	0	101	101	101	0	0	0	0;
    0	1	0	3	0	3	0	1	0	3	0	3	0	3	0	3	0	1	0	1	0	0	3	0	0	0	0	0;
    0	2	0	10	0	10	0	2	0	10	0	11	0	15	0	18	0	2	0	2	0	0	12	0	0	0	0	0;
    0	2	0	10	0	10	0	2	0	10	0	11	0	15	0	18	0	2	0	2	0	0	12	0	0	0	0	0;
    0	0	0	0	0	0	0	0	0	10	0	11	0	15	0	18	0	0	0	0	0	0	12	0	0	0	0	0];




if ia(id_b,id_a) == 1
    ia_class = IA_HEAT11();
elseif ia(id_b,id_a) == 2
    ia_class = IA_HEAT11_WATER10();
elseif ia(id_b,id_a) == 3
    ia_class = IA_HEAT11_WATER01();    
elseif ia(id_b,id_a) == 4
    ia_class = IA_HEAT11_WATER11();
elseif ia(id_b,id_a) == 5
    ia_class = IA_HEAT11_SALT10();
elseif ia(id_b,id_a) == 6
    ia_class = IA_HEAT11_SALT01();
elseif ia(id_b,id_a) == 7
    ia_class = IA_HEAT11_LAKE();
elseif ia(id_b,id_a) == 8
    ia_class = IA_HEAT11_WATER10_LAKE();
elseif ia(id_b,id_a) == 9
    ia_class = IA_HEAT11_WATER11_LAKE();
elseif ia(id_b,id_a) == 10
    ia_class = IA_HEAT11_WATER11_SNOW();
elseif ia(id_b,id_a) == 11
    ia_class = IA_HEAT11_WATER11_SNOW_XICE();
elseif ia(id_b,id_a) == 12
    ia_class = IA_HEAT11_WATER11_SNOW_LAKE();
elseif ia(id_b,id_a) == 13
    ia_class = IA_HEAT11_WATER11_LAKE_XICE();
elseif ia(id_b,id_a) == 15
    ia_class = IA_HEAT11_WATER11_RichardsEq_SNOW(); 
elseif ia(id_b,id_a) == 16
    ia_class = IA_HEAT11_WATER11_RichardsEq_LAKE(); 
elseif ia(id_b,id_a) == 18
    ia_class = IA_HEAT11_WATER11_RichardsEq_SNOW_XICE(); 
elseif ia(id_b,id_a) == 19
    ia_class = IA_HEAT11_WATER11_RichardsEq_LAKE_XICE(); 
elseif ia(id_b,id_a) == 101
    ia_class = IA_LAKE_simple_frozen_unfrozen();
else
    ia_class=0;
    disp('combination of classes not supported')
end

end