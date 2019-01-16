function stratigraphy_list = read_stratigraphies_from_file(status_info)
    %first argument: filename
    %second argiment (optional): index



    stratigraphy_list ={};
    pos_list = get_range(status_info, 'STRATIGRAPHY');
    for i=1:size(pos_list,1)
            section = status_info(pos_list(i,1):pos_list(i,2),:);
            class_name=section{2,1};
            class_handle = str2func(section{2,1});
            strat=class_handle();
            strat = initalize_from_file(strat, section);
            index = section{2,2};
            stratigraphy_list=[stratigraphy_list; {strat index}];

    end
 
end
