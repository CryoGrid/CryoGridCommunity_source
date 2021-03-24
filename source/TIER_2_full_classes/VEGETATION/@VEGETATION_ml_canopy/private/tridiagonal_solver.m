function [u] = tridiagonal_solver (a, b, c, r, n) %(a, b, c, d, n)

% Solve for U given the set of equations R * U = D, where U is a vector
% of length N, D is a vector of length N, and R is an N x N tridiagonal
% matrix defined by the vectors A, B, C each of length N. A(1) and
% C(N) are undefined and are not referenced.
%
%     |B(1) C(1) ...  ...  ...                     |
%     |A(2) B(2) C(2) ...  ...                     |
% R = |     A(3) B(3) C(3) ...                     |
%     |                    ... A(N-1) B(N-1) C(N-1)|
%     |                    ... ...    A(N)   B(N)  |
%
% The system of equations is written as:
%
%    A_i * U_i-1 + B_i * U_i + C_i * U_i+1 = D_i
%
% for i = 1 to N. The solution is found by rewriting the
% equations so that:
%
%    U_i = F_i - E_i * U_i+1

% --- Forward sweep (1 -> N) to get E and F

% % % % Matlab
% % % e(1) = c(1) / b(1);
% % % 
% % % for i = 2: 1: n-1
% % %    e(i) = c(i) / (b(i) - a(i) * e(i-1));
% % % end
% % % 
% % % f(1) = d(1) / b(1);
% % % 
% % % for i = 2: 1: n
% % %    f(i) = (d(i) - a(i) * f(i-1)) / (b(i) - a(i) * e(i-1));
% % % end
% % % 
% % % % --- Backward substitution (N -> 1) to solve for U
% % % 
% % % u(n) = f(n);
% % % 
% % % for i = n-1: -1: 1
% % %    u(i) = f(i) - e(i) * u(i+1);
% % % end


% FORTRAN
u = zeros(size(a));

bet = b(1);

u(1) = r(1) / bet;
for j = 2:n
    gam(j) = c(j-1) / bet;
    bet = b(j) - a(j)*gam(j);
    u(j) = (r(j) - a(j)*u(j-1)) / bet;
end
for j = n-1:-1:1
    u(j) = u(j) - gam(j+1)*u(j+1);
end

end


%     !
%     ! !DESCRIPTION:
%     ! Solve a tridiagonal system of equations
%     !
%     ! !USES:
%     !
%     ! !ARGUMENTS:
%     implicit none
%     integer,  intent(in)  :: n        ! Number of soil layers
%     real(r8), intent(in)  :: a(n)     ! A vector for tridiagonal solution
%     real(r8), intent(in)  :: b(n)     ! B vector for tridiagonal solution
%     real(r8), intent(in)  :: c(n)     ! C vector for tridiagonal solution
%     real(r8), intent(in)  :: r(n)     ! R vector for tridiagonal solution
%     real(r8), intent(out) :: u(n)     ! U vector for tridiagonal solution
%     !
%     ! !LOCAL VARIABLES:
%     real(r8) :: gam(n)                ! Temporary calculation
%     real(r8) :: bet                   ! Temporary calculation
%     integer  :: j                     ! Soil layer index
%     !---------------------------------------------------------------------
% 
%     ! Tridiagonal solution:
%     !
%     ! Solve for U given the set of equations F x U = R, where U is a vector
%     ! of length N, R is a vector of length N, and F is an N x N tridiagonal
%     ! matrix defined by the vectors A, B, C (each of length N). A(1) and
%     ! C(N) are undefined and are not referenced by the subroutine.
%     !
%     !    | b(1) c(1)   0  ...                      |   | u(1)   |   | r(1)   |
%     !    | a(2) b(2) c(2) ...                      |   | u(2)   |   | r(2)   |
%     !    |                ...                      | x | ...    | = | ...    |
%     !    |                ... a(n-1) b(n-1) c(n-1) |   | u(n-1) |   | r(n-1) |
%     !    |                ...   0    a(n)   b(n)   |   | u(n)   |   | r(n)   |
%     !


