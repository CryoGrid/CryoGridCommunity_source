function forcing_ini_inputs_check( ini_struct )
% This function verifies the validity of the potential variables' inputs
% extracted from the resulting structure of the ini parser.

ini_forcing_name = fieldnames(ini_struct);
ini_forcing_variables = fields(ini_struct.(ini_forcing_name{1,1}));

%Check the data type, extension, range, etc. of the different inputs from
%potentially existing variables

if ~isempty(find(strcmp(ini_forcing_variables, 'filename') == 1))
    i1 = find(strcmp(ini_forcing_variables, 'filename') == 1);
    validateattributes(ini_struct.(ini_forcing_name{1}).(ini_forcing_variables{i1,1}), {'char'},{'scalartext'});
    if ~contains(ini_struct.(ini_forcing_name{1}).(ini_forcing_variables{i1,1}),'.mat')
        error('The expected extension of the forcing dataset is .mat');
    end
end

if ~isempty(find(strcmp(ini_forcing_variables, 'start_time') == 1))
    i2 = find(strcmp(ini_forcing_variables, 'start_time') == 1);
    is_valid_date_ddmmyyyy(ini_struct.(ini_forcing_name{1}).(ini_forcing_variables{i2,1}));
end

if ~isempty(find(strcmp(ini_forcing_variables, 'end_time') == 1))
    i3 = find(strcmp(ini_forcing_variables, 'end_time') == 1);
    is_valid_date_ddmmyyyy(ini_struct.(ini_forcing_name{1}).(ini_forcing_variables{i3,1}));
end

if ~isempty(find(strcmp(ini_forcing_variables, 'rain_fraction') == 1))
    i4 = find(strcmp(ini_forcing_variables, 'rain_fraction') == 1);
    is_boolean(ini_struct.(ini_forcing_name{1}).(ini_forcing_variables{i4,1}))
end

if ~isempty(find(strcmp(ini_forcing_variables, 'snow_fraction') == 1))
    i5 = find(strcmp(ini_forcing_variables, 'snow_fraction') == 1);
    is_boolean(ini_struct.(ini_forcing_name{1}).(ini_forcing_variables{i5,1}))
end

if ~isempty(find(strcmp(ini_forcing_variables, 'latitude') == 1))
    i6 = find(strcmp(ini_forcing_variables, 'latitude') == 1);
    validateattributes(ini_struct.(ini_forcing_name{1}).(ini_forcing_variables{i6,1}), {'numeric'},{'>=',-90,'<=',90});
end

if ~isempty(find(strcmp(ini_forcing_variables, 'longitude') == 1))
    i7 = find(strcmp(ini_forcing_variables, 'longitude') == 1);
    validateattributes(ini_struct.(ini_forcing_name{1}).(ini_forcing_variables{i7,1}), {'numeric'},{'>=',-180,'<=',180});
end

if ~isempty(find(strcmp(ini_forcing_variables, 'altitude') == 1))
    i8 = find(strcmp(ini_forcing_variables, 'altitude') == 1);
    validateattributes(ini_struct.(ini_forcing_name{1}).(ini_forcing_variables{i8,1}), {'numeric'},{'>=',0,'<',8900});
end

if ~isempty(find(strcmp(ini_forcing_variables, 'domain_depth') == 1 ))
    i9 = find(strcmp(ini_forcing_variables, 'domain_depth') == 1 );
    validateattributes(ini_struct.(ini_forcing_name{1}).(ini_forcing_variables{i9,1}), {'numeric'},{'positive'});
end

if ~isempty(find(strcmp(ini_forcing_variables, 'heatFlux_lb') == 1))
    i10 = find(strcmp(ini_forcing_variables, 'heatFlux_lb') == 1);
    validateattributes(ini_struct.(ini_forcing_name{1}).(ini_forcing_variables{i10,1}), {'single','double'},{'>',0,'<',0.11});
end

if ~isempty(find(strcmp(ini_forcing_variables, 'airT_height') == 1))
    i11 = find(strcmp(ini_forcing_variables, 'airT_height') == 1);
    validateattributes(ini_struct.(ini_forcing_name{1}).(ini_forcing_variables{i11,1}), {'numeric'},{'>=',0});
end

%i12, etc. if more variables are necessary to initialize the different
%forcing classes

end

