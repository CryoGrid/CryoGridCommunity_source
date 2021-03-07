
num_variables = 19;
num_ranks = 100;

res=[];

for my_rank=1:100;

core_list = ceil([0:num_variables./num_ranks:num_variables])';
%core_list = max(1, core_list_raw);
core_list = core_list(2:end,1);
variable_start_index = core_list(my_rank,1);
permission_list = zeros(num_variables,1);
dummy = [0; core_list];
if dummy(my_rank+1,1) ~= dummy(my_rank,1)
    permission_list(variable_start_index,1)=1;
    if my_rank == num_ranks
        permission_list(variable_start_index:num_variables,1) = 1;
        permission_list(1:core_list(1,1)-1,1) = 1;
    else
        permission_list(variable_start_index:core_list(my_rank+1,1)-1,1) = 1;
    end
end

res=[res permission_list];
end