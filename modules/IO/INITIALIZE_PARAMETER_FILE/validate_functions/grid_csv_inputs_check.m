function grid_csv_inputs_check ( csv_table )
% This function verifies the validity of the  inputs
% extracted from the resulting table of the csv parser.

% Check whether upper boundaries in the grid table are integer values
is_integer_value(csv_table.upper);

% Check whether spacings in the grid table are numeric values
validateattributes(csv_table.spacing, {'numeric'},{'positive'});

% Check whether lower boundaries in the grid table are integer values and
% correspond to the next upper boundaries in the grid table
is_integer_value(csv_table.lower);
for i = 1:(length(csv_table.lower)-1)
    if csv_table.lower(i) ~= csv_table.upper
        error('The lower boundary in grid.lower should be equal to the next upper boundary in grid.upper')
    end
end

end

