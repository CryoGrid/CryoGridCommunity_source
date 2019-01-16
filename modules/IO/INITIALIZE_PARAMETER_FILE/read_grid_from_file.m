function grid = read_grid_from_file(varargin)
    %first argument: filename
    %second argiment (optional): index
    status_info = varargin{1};
    if nargin == 1
        index=1;
    else
        index = varargin{2};
    end

    
    pos_list = get_range(status_info, 'GRID');
    for i=1:size(pos_list,1)
        if status_info{pos_list(i,1)+1,2} == index 
            section = status_info(pos_list(i,1):pos_list(i,2),:);
            class_name=section{2,1};
            class_handle = str2func(section{2,1});
            grid=class_handle();
            grid = initalize_from_file(grid, section);

        end
    end
 
end
