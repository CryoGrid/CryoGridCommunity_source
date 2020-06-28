function is_integer_value( input_value )
% This function checks whether the input is an integer.

if isa(input_value,'char') == 1 | any(mod(input_value,1)~= 0) == 1
    error('Integer data type is expected here')
end
        
end

