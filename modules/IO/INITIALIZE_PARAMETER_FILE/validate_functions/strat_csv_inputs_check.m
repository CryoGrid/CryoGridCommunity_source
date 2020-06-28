function strat_csv_inputs_check ( csv_table, class_name, available_classes_list )
% This function verifies the validity of the  inputs
% extracted from the resulting table of the csv parser. As the variables
% differ greatly from one STRAT class to the other, the inputs being
% checked will be dependent on the name of the class.

if strcmp(class_name, 'STRAT_linear') == 1
    
    if isnumeric(csv_table.T) == 0
        error('Positive, negative (or null) numeric values are expected')
    end
    
    validateattributes(csv_table.depth, {'numeric'},{'>=',0});

elseif strcmp(class_name, 'STRAT_layers') == 1
    
    validateattributes(csv_table.depth, {'numeric'},{'>=',0});
    
    validateattributes(csv_table.mineral, {'numeric'},{'>=',0,'<=',1}); 
    validateattributes(csv_table.organic, {'numeric'},{'>=',0,'<=',1});
    validateattributes(csv_table.waterIce, {'numeric'},{'>=',0,'<=',1});
    % Not certain about the range and unit of Xice, to verify !
    validateattributes(csv_table.Xice, {'numeric'},{'>=',0,'<=',1});  
    if any((csv_table.mineral + csv_table.organic + csv_table.waterIce + csv_table.Xice) ~= 1) == 1
        error('The total sum of the mineral, organic, waterIce and Xice fractions in a soil layer should equal 1')
    end
    
    % Lines should be added to check the data type and range of the
    % soil_type_code variable, but not used for now in the main routine.
    % The unit and range of the field_capacity variable should be verified
    % before adding the corresponding validate functions.
    
    validateattributes(csv_table.saltConc, {'numeric'}, {'>=',0});
    
    validateattributes(csv_table.deltaT, {'numeric'}, {'>=',0})

elseif strcmp(class_name, 'STRAT_classes') == 1
    
    validateattributes(csv_table.depth, {'numeric'}, {'>=',0});
    
    for i = 1:length(csv_table.class_name)
        class_exists(csv_table.class_name{i}, available_classes_list)
    end
    for i = 1:(length(csv_table.snow_class)-1)
        class_exists(csv_table.snow_class{i}, available_classes_list)
    end
    
    is_positive_integer_value(csv_table.groundindex);
    for i = 1:(length(csv_table.snowindex)-1)
        is_positive_integer_value(csv_table.snowindex(i));
    end
end

