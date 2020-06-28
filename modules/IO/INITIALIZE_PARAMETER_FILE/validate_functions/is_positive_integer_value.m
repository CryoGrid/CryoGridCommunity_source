function is_positive_integer_value( input_value )
% This function checks whether the input is an integer and is strictly positive.

if isa(input_value,'char') == 1 | any(mod(input_value,1)~= 0) == 1 | any(input_value <= 0) ==1
    error('Positive integer data type is expected here')
end
        
end

