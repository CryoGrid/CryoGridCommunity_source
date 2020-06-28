function class_list = read_classes_from_file(status_info)
    %first argument: filename
    %second argiment (optional): index

    
    
    
    
    

    class_list ={};
    pos_list = get_range(status_info, 'CLASS');
    for i=1:size(pos_list,1)
            section = status_info(pos_list(i,1):pos_list(i,2),:);
            class_name=section{2,1};
            class_handle = str2func(section{2,1});
            class = class_handle();
            class = provide_variables(class);
            class.PARA = initialize_from_file(class, class.PARA, section);
            index = section{2,2};
            class_list = [class_list; {class index}];
    end
 
end
