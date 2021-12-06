function x = TDMAsolver(a,b,c,d)

%a, b, c are the column vectors for the compressed tridiagonal matrix, d is the right vector
n = length(b); % n is the number of rows
x = zeros(n,1);
% Modify the first-row coefficients
c(1) = c(1) / b(1);    % Division by zero risk.
d(1) = d(1) / b(1);    % Division by zero would imply a singular matrix.

for i = 2:n-1
    temp = b(i) - a(i) * c(i-1);
    c(i) = c(i) / temp;
    d(i) = (d(i) - a(i) * d(i-1)) / temp;
end

d(n) = (d(n) - a(n) * d(n-1))/( b(n) - a(n) * c(n-1));

% Now back substitute.
x(n) = d(n);
for i = n-1:-1:1
    x(i) = d(i) - c(i) * x(i + 1);
end

end
