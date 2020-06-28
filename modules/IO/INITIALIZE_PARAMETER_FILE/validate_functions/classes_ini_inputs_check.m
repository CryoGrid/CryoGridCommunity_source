function classes_ini_inputs_check( ini_struct )
% This function verifies the validity of the potential variables' inputs
% extracted from the resulting structure of the ini parser.

ini_classes_names = fieldnames(ini_struct);

%Check the data type, extension, range, etc. of the different inputs from
%potentially existing variables

for i = 1:length(ini_classes_names)
	ini_classes_variables = fields(ini_struct.(ini_classes_names{i,1}));
	
	%This check needs some improvements to handle arrays
	if ~isempty(find(strcmp(ini_classes_variables, 'albedo') == 1))
		i1 = find(strcmp(ini_classes_variables, 'albedo') == 1);
		validateattributes(ini_struct.(ini_classes_names{i,1}).(ini_classes_variables{i1,1}), {'numeric'},{'>=',0,'<=',1});
	end

	if ~isempty(find(strcmp(ini_classes_variables, 'epsilon') == 1))
		i2 = find(strcmp(ini_classes_variables, 'epsilon') == 1);
		validateattributes(ini_struct.(ini_classes_names{i,1}).(ini_classes_variables{i2,1}), {'numeric'},{'>=',0,'<',1});
	end

	if ~isempty(find(strcmp(ini_classes_variables, 'rs') == 1))
		i3 = find(strcmp(ini_classes_variables, 'rs') == 1);
		validateattributes(ini_struct.(ini_classes_names{i,1}).(ini_classes_variables{i3,1}), {'numeric'},{'>=',0});
	end
	
	if ~isempty(find(strcmp(ini_classes_variables, 'z0') == 1))
		i4 = find(strcmp(ini_classes_variables, 'z0') == 1);
		validateattributes(ini_struct.(ini_classes_names{i,1}).(ini_classes_variables{i4,1}), {'numeric'},{'>=',0,'<',2});
	end

	if ~isempty(find(strcmp(ini_classes_variables, 'dt_max') == 1))
		i5 = find(strcmp(ini_classes_variables, 'dt_max') == 1);
		validateattributes(ini_struct.(ini_classes_names{i,1}).(ini_classes_variables{i5,1}), {'numeric'},{'positive'});
	end
	
	if ~isempty(find(strcmp(ini_classes_variables, 'dE_max') == 1))
		i6 = find(strcmp(ini_classes_variables, 'dE_max') == 1);
		validateattributes(ini_struct.(ini_classes_names{i,1}).(ini_classes_variables{i6,1}), {'numeric'},{'positive'});
	end
	
	if ~isempty(find(strcmp(ini_classes_variables, 'density') == 1))
		i7 = find(strcmp(ini_classes_variables, 'density') == 1);
		validateattributes(ini_struct.(ini_classes_names{i,1}).(ini_classes_variables{i7,1}), {'numeric'},{'>=',0,'<',1000});
	end
	
	if ~isempty(find(strcmp(ini_classes_variables, 'field_capacity') == 1))
		i8 = find(strcmp(ini_classes_variables, 'field_capacity') == 1);
		validateattributes(ini_struct.(ini_classes_names{i,1}).(ini_classes_variables{i8,1}), {'numeric'},{'>=',0,'<=',1});
	end
	
	if ~isempty(find(strcmp(ini_classes_variables, 'swe_per_cell') == 1))
		i9 = find(strcmp(ini_classes_variables, 'swe_per_cell') == 1);
		validateattributes(ini_struct.(ini_classes_names{i,1}).(ini_classes_variables{i9,1}), {'numeric'},{'>=',0,'<=',1});
	end
	
	% The range of validity might need to be detailed
	if ~isempty(find(strcmp(ini_classes_variables, 'tortuosity') == 1))
		i10 = find(strcmp(ini_classes_variables, 'tortuosity') == 1);
		validateattributes(ini_struct.(ini_classes_names{i,1}).(ini_classes_variables{i10,1}), {'numeric'},{'positive'});
    end
    
    if ~isempty(find(strcmp(ini_classes_variables, 'albedo_max') == 1))
		i11 = find(strcmp(ini_classes_variables, 'albedo_max') == 1);
		validateattributes(ini_struct.(ini_classes_names{i,1}).(ini_classes_variables{i11,1}), {'numeric'},{'>=',0,'<=',1});
    end
    
     if ~isempty(find(strcmp(ini_classes_variables, 'albedo_min') == 1))
		i12 = find(strcmp(ini_classes_variables, 'albedo_min') == 1);
		validateattributes(ini_struct.(ini_classes_names{i,1}).(ini_classes_variables{i12,1}), {'numeric'},{'>=',0,'<=',1});
     end
    
     if ~isempty(find(strcmp(ini_classes_variables, 'albedo_max') == 1)) & ~isempty(find(strcmp(ini_classes_variables, 'albedo_min') == 1))
        i11 = find(strcmp(ini_classes_variables, 'albedo_max') == 1);
        i12 = find(strcmp(ini_classes_variables, 'albedo_min') == 1);
        if ini_struct.(ini_classes_names{i,1}).(ini_classes_variables{i12,1}) > ini_struct.(ini_classes_names{i,1}).(ini_classes_variables{i11,1})
            error('The variable albedo_min must be given a value inferior or equal to albedo_max and vice-versa')
        end
     end
end

end
	
	