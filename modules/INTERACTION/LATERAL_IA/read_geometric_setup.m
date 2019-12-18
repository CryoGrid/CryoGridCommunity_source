function lateral = read_geometric_setup(lateral,result_path,run_number)

    load([result_path '\' run_number '\geometric_setup.mat']);
    
    lateral.PARA.relative_elevation = geometry.relative_elevation;
    lateral.PARA.area = geometry.area;
    lateral.PARA.distance = geometry.distance;
    lateral.PARA.contact_length = geometry.contact_length;
    lateral.PARA.exposure = geometry.exposure;
    
end