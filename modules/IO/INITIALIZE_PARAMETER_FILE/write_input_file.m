function classes= write_input_file(varargin)

addpath('modules')

disp("Number of input arguments: " + nargin)
 %   celldisp(varargin)

    classes={};

for i=1:size(varargin,2)
    %classes = 
    dummy=varargin{i};
    classes=[classes; dummy()];
    
end

% [a,b,c]=xlsread('test.xls');
% 
% A={'run_1_test'; ''; 'albedo'; 'z0'};
% xlswrite('test.xls', A);