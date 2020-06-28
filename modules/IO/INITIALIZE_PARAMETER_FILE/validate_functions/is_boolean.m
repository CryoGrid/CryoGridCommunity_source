function is_boolean( input_value )
%This function checks whether the input is of type boolean, meaning if it
%has a value of either 0 or 1.

if input_value ~= 0 & input_value ~= 1
    error('Values of 1 or 0 are expected')
end

end

