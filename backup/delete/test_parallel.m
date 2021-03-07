load_index = [1; 1; 1.8];
number_of_cores = 5;
old_end_time = 0;

for my_core = 1:number_of_cores
    %my_core = 12;
active = [0;0;0];
number_of_years = [];

number_of_slices = 46;
start_year= [];

for i=1:size(run_info.PARA.tile_preproc_class_index,1)
    tile = run_info.PPROVIDER.CLASSES.(run_info.PARA.tile_preproc_class{i,1}){run_info.PARA.tile_preproc_class_index(i,1),1};
    forcing = run_info.PPROVIDER.CLASSES.(tile.PARA.forcing_class){tile.PARA.forcing_class_index,1};
    number_of_years = [number_of_years; forcing.PARA.end_time(1)-forcing.PARA.start_time(1)+1];
    start_year= [start_year; forcing.PARA.start_time(1)];
end

load_per_core = sum(number_of_years .* number_of_slices .* load_index) ./ number_of_cores;
number_of_slices_per_tile = number_of_years .* number_of_slices;
load_per_tile = number_of_years .* number_of_slices .* load_index;
load_per_tile_acc = cumsum(load_per_tile);

result =[];
tile_number = 1;
core = 1;
slice_count = 0;
end_last_slice = 0;
load_count = 0;
slice_count_prev = 1;


load_start = load_per_core .* (my_core-1)+1;
load_end = load_per_core .* my_core;

tile_num = 1;
slice_count = 0; 
while tile_num<=3 
    if load_start > load_per_tile_acc(tile_num) 
       tile_num=tile_num+1;
    else
        break
    end

end
tile_num_start = tile_num;

tile_num = 1;
slice_count = 0; 
while tile_num<=3 
    if load_end > load_per_tile_acc(tile_num) 
       tile_num=tile_num+1;
    else
        break
    end

end
tile_num_end = tile_num;
active(tile_num_start:tile_num_end) = active(tile_num_start:tile_num_end)+1;


load_per_tile_acc = [0; load_per_tile_acc];
slice_count_start = round((load_start - load_per_tile_acc(tile_num_start)) ./ load_index(tile_num_start));
if slice_count_start > 1
%     slice_count_start
    new_start_time = datenum(start_year(tile_num_start,1) + floor((slice_count_start-1)./46),1,1) + (mod(slice_count_start-1, 46)).*8;
    disp(['start ' datestr(new_start_time)])
end
slice_count_end =  round((load_end - load_per_tile_acc(tile_num_end)) ./ load_index(tile_num_end))+2; %add the 2 for overlap
if slice_count_end < number_of_slices_per_tile(tile_num_end)
    new_end_time = datenum(start_year(tile_num_end,1) + floor((slice_count_end-1)./46),1,1) + (mod(slice_count_end-1, 46)).*8;
    datestr(new_end_time)
    disp(old_end_time - new_start_time)
end
old_end_time = new_end_time;
end

% load_end >= load_per_tile_acc(tile_num)

% 
% number_of_slices_acc = cumsum(number_of_slices_per_tile);
% tile_num=1;
% while result(my_core,3) > number_of_slices_acc(tile_num)
%     tile_num = tile_num+1;
% end
% active(tile_num)=1;
% if tile_num==1
%     new_start_time = datenum(start_year(tile_num,1) + floor((result(tile_num,3))./46),1,1) + (mod(result(tile_num,3),46)-1).*8;
% else
%     new_start_time = datenum(start_year(tile_num,1) + floor((result(tile_num,3)-number_of_slices_acc(tile_num-1,1))./46),1,1) + mod(result(tile_num,3),46).*8-1;
% end
% datestr(new_start_time)
% 
% while result(my_core,4) > number_of_slices_acc(tile_num)
%     
%     tile_num = tile_num+1;
%     active(tile_num)=1;
% end
% 
% 
% 
% % while core <= number_of_cores %tile_number<=3 
% %         load_count = load_count + load_per_core;
% %         if load_count < load_per_tile_acc(tile_number)
% %             slice_count = slice_count + load_per_core./load_index(tile_number);
% %             result=[result; [core tile_number round(slice_count)]];
% %         else
% %             while tile_number<3
% %                 remaining_load = load_count - load_per_tile_acc(tile_number);
% %                 if  remaining_load >0
% %                     result=[result; [core tile_number number_of_slices_per_tile(tile_number)]];
% %                 else
% %                     break
% %                 end
% %                 tile_number = tile_number+1;
% %                 slice_count =  remaining_load./load_index(tile_number);
% %                 if round(slice_count) <= number_of_slices_per_tile(tile_number) 
% %                     result=[result; [core tile_number round(slice_count)]];
% %                 end
% %             end
% %         end
% %         
% %         core = core +1;
% % end

% result2=[];
% ti=0;
% for i=1:size(result,1)
%     
%     if ti ~= result(i,2)
%         result2 = [result2; [start_year(result(i,2)) + floor(result(i,3)./46) mod(result(i,3),46).*8 + 1 ...
%         datenum(start_year(result(i,2)) + floor(result(i,3)./46),1,1) + mod(result(i,3),46).*8-1  ...
%         datenum(start_year(result(i,2)),1,1)]];
%         ti = result(i,2);
%     else
%         result2 = [result2; [start_year(result(i,2)) + floor(result(i,3)./46) mod(result(i,3),46).*8 + 1 ...
%         datenum(start_year(result(i,2)) + floor(result(i,3)./46),1,1) + mod(result(i,3),46).*8-1  ...
%         datenum(start_year(result(i-1,2)) + floor(result(i-1,3)./46),1,1) + mod(result(i-1,3),46).*8]];
%     end
% end
% % result2 = [result2 [datenum(start_year(1,1), 1, 1); result2(1:end-1,3)]];
% % result2(:,3) = result2(:,3)-1;
% 
% for i=1:size(result,1)
%    if my_core == result(i,1)
%        active(result(i,2),1) = 1;
%        new_start_time = datestr(result2(i,4))
%        new_end_time = datestr(result2(i,3))
%    end
%     
% end
    
