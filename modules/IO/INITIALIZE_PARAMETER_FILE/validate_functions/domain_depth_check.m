function  domain_depth_check ( grid, forcing )
% This function verifies whether the value of the domain_depth variable set in forcing matches one of the grid point.
% If it is not the case, an error message is displayed.

if isempty(find(grid.GRID == forcing.PARA.domain_depth)) == 1
	error ('domain_depth in forcing is expected to match one of the values of the grid. Please attribute a valid value (in m) to the domain_depth variable in FORCING.ini')
end

end