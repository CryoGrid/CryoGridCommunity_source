%Mathtools%-----------------------------------------------------------------------
function [r1, r2] = quadratic(a, b, c) % (a, b, c, r1, r2)

% DESCRIPTION:
% Solve a quadratic equation for its two roots

%ARGUMENTS:
%   a, b, c       % Terms for quadratic equation
%   r1, r2       % Roots of quadratic equation

%LOCAL VARIABLES:
%   q                        % Temporary term for quadratic solution


% Whatever values given as a, b, c, right before the function is executed.
%  a = 1; % dummy values
%  b = 1; % dummy values
%  c = 1; % dummy values


if (a == 0)
    fprintf("Quadratic solution error: a = ", a); %?
    %     call endrun()
end

if (b >= 0 )
    q = -0.5  .* (b + sqrt(b.*b - 4 .*a.*c));
else
    q = -0.5  .* (b - sqrt(b.*b - 4 .*a.*c));
end

r1 = q ./ a;
if (q ~= 0)  
    r2 = c ./ q;
else
    r2 = 1.e36;
end
end
