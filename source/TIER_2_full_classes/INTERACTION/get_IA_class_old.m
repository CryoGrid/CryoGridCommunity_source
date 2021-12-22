%========================================================================
% CryoGrid get_IA_class function
% identifies and initializes the correct interaction class between pairs of
% GROUND classes, returns 0 if classes are not compatible
% INPUT: GROUND class nacmes as strings
% NOTE: be very careful when updating the interaction matrix "ia"! Keep the
% last working version in the code (commented out)
% S. Westermann, October 2020 - discontinued Dec 2021!!!
%========================================================================


function ia_class = get_IA_class_old(above_class, below_class)

%---------------
if strcmp(above_class(end-8:end), 'FLIP_FLOP')
    above_class = above_class(1:end-10);
end

if strcmp(below_class(end-8:end), 'FLIP_FLOP')
    below_class = below_class(1:end-10);
end

%quick fixes
%---------
if strcmp(above_class, 'GROUND_freezeC_bucketW_Xice_seb_snow_BGC')
    above_class = 'GROUND_freezeC_bucketW_Xice_seb_snow';
end
if strcmp(below_class, 'GROUND_freezeC_bucketW_Xice_seb_snow_BGC')
    below_class = 'GROUND_freezeC_bucketW_Xice_seb_snow';
end


if strcmp(above_class, 'GROUND_freezeC_RichardsEqW_seb_pressure')
    above_class = 'GROUND_freezeC_RichardsEqW_seb';
end
if strcmp(below_class, 'GROUND_freezeC_RichardsEqW_seb_pressure')
    below_class = 'GROUND_freezeC_RichardsEqW_seb';
end

if strcmp(above_class, 'GROUND_freezeC_RichardsEqW_seb_pressure_snow')
    above_class = 'GROUND_freezeC_RichardsEqW_seb_snow';
end
if strcmp(below_class, 'GROUND_freezeC_RichardsEqW_seb_pressure_snow')
    below_class = 'GROUND_freezeC_RichardsEqW_seb_snow';
end

%----------
if strcmp(above_class, 'GROUND_freezeC_ubT')
    above_class = 'GROUND_freezeC_seb';
end
if strcmp(below_class, 'GROUND_freezeC_ubT')
    below_class = 'GROUND_freezeC_seb';
end

if strcmp(above_class, 'GROUND_freezeC_ubT_snow')
    above_class = 'GROUND_freezeC_seb_snow';
end
if strcmp(below_class, 'GROUND_freezeC_ubT_snow')
    below_class = 'GROUND_freezeC_seb_snow';
end

%----------
if strcmp(above_class, 'GROUND_freeW_ubT')
    above_class = 'GROUND_freeW_seb';
end
if strcmp(below_class, 'GROUND_freeW_ubT')
    below_class = 'GROUND_freeW_seb';
end

if strcmp(above_class, 'GROUND_freeW_ubT_snow')
    above_class = 'GROUND_freeW_seb_snow';
end
if strcmp(below_class, 'GROUND_freeW_ubT_snow')
    below_class = 'GROUND_freeW_seb_snow';
end
%------

if strcmp(above_class, 'GROUND_fcSimple_salt_ubT')
    above_class = 'GROUND_fcSimple_salt_seb';
end
if strcmp(below_class, 'GROUND_fcSimple_salt_ubT')
    below_class = 'GROUND_fcSimple_salt_seb';
end

if strcmp(above_class, 'GROUND_fcSimple_salt_ubT_snow')
    above_class = 'GROUND_fcSimple_salt_seb_snow';
end
if strcmp(below_class, 'GROUND_fcSimple_salt_ubT_snow')
    below_class = 'GROUND_fcSimple_salt_seb_snow';
end
%------

if strcmp(above_class, 'SNOW_simple_ubT')
    above_class = 'SNOW_simple_seb';
end
if strcmp(below_class, 'SNOW_simple_ubT')
    below_class = 'SNOW_simple_seb';
end


%------------------

if strcmp(above_class, 'GROUND_store_flip_flop_singleClass') || strcmp(above_class, 'GROUND_store_flip_flop_singleClass_BGC')
    above_class = 'GROUND_freeW_seb';
end


if strcmp(above_class, 'GLACIER_freeW_seb_snow')
    above_class = 'GROUND_freeW_seb_snow';
end
if strcmp(below_class, 'GLACIER_freeW_seb_snow')
    below_class = 'GROUND_freeW_seb_snow';
end

%set dummy, so that it doesn't crash
if strcmp(above_class, 'GROUND_multi_tile2') || strcmp(above_class, 'GROUND_multi_tile2_Cas')
    above_class2 = above_class;
    above_class = 'GROUND_freeW_seb';
else
    above_class2 = 0;
end



%list of all CryoGrid GROUND classes
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
    {'GROUND_freezeC_RichardsEqW_seb_vegetation'};
    {'GROUND_freezeC_RichardsEqW_seb_vegetation_snow'};
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
    {'SNOW_crocus_bucketW_seb_vegetation'};
    {'SNOW_crocus2_bucketW_seb'}];



for i = 1:length(classes)
    if strcmp(below_class,classes{i})
        id_a = i;
    end
    if strcmp(above_class,classes{i})
        id_b = i;
    end
end




%interaction matrix
% ia=[1	0	3	3	3	3	1	0	3	3	3	3	0	0	0	0	6	6	0	0	0	0	0	0	0	0	0	0;
%     1	1	3	3	3	3	1	1	3	3	3	3	0	0	0	0	6	6	0	0	0	0	0	0	0	0	0	0;
%     2	0	4	4	4	4	2	0	4	4	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
%     2	0	4	0	4	0	2	0	4	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
%     2	0	4	4	0	0	2	0	4	4	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
%     2	0	4	0	0	0	2	0	4	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
%     1	0	3	3	3	3	1	0	3	3	3	3	0	0	0	0	6	6	0	0	0	0	0	0	0	0	0	0;
%     1	1	3	3	3	3	1	1	3	3	3	3	0	0	0	0	6	6	0	0	0	0	0	0	0	0	0	0;
%     2	0	4	4	4	4	2	0	4	4	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
%     2	0	4	0	4	0	2	0	4	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
%     2	0	0	0	0	0	2	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
%     2	0	0	0	0	0	2	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
%     2	2	0	0	0	0	2	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
%     2	2	0	0	0	0	2	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
%     2	2	0	0	0	0	2	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
%     2	2	0	0	0	0	2	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
%     5	0	0	0	0	0	5	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
%     5	0	0	0	0	0	5	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
%     1	1	0	0	0	0	1	1	0	0	0	0	0	0	0	0	0	0	0	0	101	0	0	0	0	0	0	0;
%     1	1	0	0	0	0	1	1	0	0	0	0	0	0	0	0	0	0	0	0	101	0	0	0	0	0	0	0;
%     7	7	0	0	0	0	7	7	0	0	0	0	0	0	0	0	0	0	101	101	0	0	0	0	0	0	0	0;
%     8	8	9	9	9	9	8	8	9	9	13	13	16	16	19	19	0	0	0	0	0	101	101	101	0	0	0	0;
%     8	8	9	9	9	9	8	8	9	9	13	13	16	16	19	19	0	0	0	0	0	101	101	101	0	0	0	0;
%     8	8	9	9	9	9	8	8	9	9	13	13	16	16	19	19	0	0	0	0	0	101	101	101	0	0	0	0;
%     0	1	0	3	0	3	0	1	0	3	0	3	0	3	0	3	0	1	0	1	0	0	3	0	0	0	0	0;
%     0	2	0	10	0	10	0	2	0	10	0	11	0	15	0	18	0	2	0	2	0	0	12	0	0	0	0	0;
%     0	2	0	10	0	10	0	2	0	10	0	11	0	15	0	18	0	2	0	2	0	0	12	0	0	0	0	0;
%     0	0	0	0	0	0	0	0	0	10	0	11	0	15	0	18	0	0	0	0	0	0	12	0	0	0	0	0];


ia=[1	0	3	3	3	3	1	0	3	3	3	3	0	0	0	0	0	0	6	6	0	0	0	0	0	0	0	0	0	0	0;
    1	1	3	3	3	3	1	1	3	3	3	3	0	0	0	0	0	0	6	6	0	0	0	0	0	0	0	0	0	0	0;
    2	0	4	4	4	4	2	0	4	4	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
    2	0	4	0	4	0	2	0	4	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
    2	0	4	4	0	0	2	0	4	4	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
    2	0	4	0	0	0	2	0	4	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
    1	0	3	3	3	3	1	0	3	3	3	3	0	0	0	0	0	0	6	6	0	0	0	0	0	0	0	0	0	0	0;
    1	1	3	3	3	3	1	1	3	3	3	3	0	0	0	0	0	0	6	6	0	0	0	0	0	0	0	0	0	0	0;
    2	0	4	4	4	4	2	0	4	4	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
    2	0	4	0	4	0	2	0	4	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
    2	0	0	0	0	0	2	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
    2	0	0	0	0	0	2	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
    2	2	0	0	0	0	2	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
    2	2	0	0	0	0	2	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
    2	2	0	0	0	0	2	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
    2	2	0	0	0	0	2	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
    2	2	0	0	0	0	2	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
    2	2	0	0	0	0	2	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
    5	0	0	0	0	0	5	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
    5	0	0	0	0	0	5	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
    1	1	0	0	0	0	1	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	101	0	0	0	0	0	0	0	0;
    1	1	0	0	0	0	1	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	101	0	0	0	0	0	0	0	0;
    7	7	0	0	0	0	7	7	0	0	0	0	0	0	0	0	0	0	0	0	101	101	0	0	0	0	0	0	0	0	0;
    8	8	9	9	9	9	8	8	9	9	13	13	16	16	16	16	19	19	0	0	0	0	0	101	101	101	0	0	0	0	0;
    8	8	9	9	9	9	8	8	9	9	13	13	16	16	16	16	19	19	0	0	0	0	0	101	101	101	0	0	0	0	0;
    8	8	9	9	9	9	8	8	9	9	13	13	16	16	16	16	19	19	0	0	0	0	0	101	101	101	0	0	0	0	0;
    0	1	0	3	0	3	0	1	0	3	0	3	0	3	0	3	0	3	0	1	0	1	0	0	3	0	0	0	0	0	0;
    0	2	0	10	0	10	0	2	0	10	0	11	0	15	0	15	0	18	0	2	0	2	0	0	12	0	0	0	0	0	0;
    0	2	0	10	0	10	0	2	0	10	0	11	0	15	0	15	0	18	0	2	0	2	0	0	12	0	0	0	0	0	0;
    0	2	0	10	0	10	0	2	0	10	0	11	0	15	0	15	0	18	0	2	0	2	0	0	12	0	0	0	0	0	0;
    0	0	0	0	0	0	0	0	0	10	0	11	0	15	0	15	0	18	0	0	0	0	0	0	12	0	0	0	0	0	0];



%all avilable INTERACTION (IA) classes
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
    %disp('combination of classes not supported')
end


if strcmp(above_class2, 'GROUND_multi_tile2') || strcmp(above_class2, 'GROUND_multi_tile2_Cas')
    ia_class = IA_MULTI_TILE10();
end
% if  strcmp(below_class,'READ_OUT_BGC') && strcmp(above_class,'READ_OUT_BGC')
%     ia_class = IA_DO_NOTHING();
% end

end