%==========================================================================
% Generate class inheritance diagrams from main CryoGrid classes Tier 2&3
% T. Ingemann-Nielsen, November 2020
%==========================================================================

% Folders to include in automatic UML chart generation
% column 1 is tag to include in beginning of UML filename
% column 2 is relative path to the classes to include (can be to a single class file, or contain wildcards)
% This definition will make one hierachy diagram for each class found in
% the TIER_2 and TIER_3 folders (with exceptions defined below).
folder_list = {'TIER2' '..\modules\TIER_2_full_classes\**\*.m';
               'TIER3' '..\modules\TIER_3_snow\**\*.m';
               '' ''};      % the last line is needed to ensure code runs also if only one folder/file is specified
           
% any class for which the path starts with one of these will be excluded
exclude_folders = {'..\modules\TIER_2_full_classes\BGC'
                   '..\modules\TIER_2_full_classes\INTERACTION'
                   };

% The use of relative paths in matlab is not well supported.
% Carefully test, if adding additional paths above.

process_Tier2_Tier3 = true;
open_svg = false;   % Set to true to open each svg file in the browser window


% MAIN CODE BELOW

try
    m2uml.create_PlantUML_script();
catch
    disp('This script requires the toolbox m2uml.')
    disp('Please download and install from https://se.mathworks.com/matlabcentral/fileexchange/59722-m2uml')
    error('Toolbox m2uml missing!')
end

if contains(pwd, ' ')
    disp(['Current working directory: ' pwd]);
    error('The m2uml toolbox does not allow spaces in paths. Please rename the path to CryoGrid');
end

addpath(genpath('..\modules'));


if process_Tier2_Tier3
    % Generate class list according to folder specified
    class_list = {};
    id = 1;

    for m = 1:length(folder_list)
        % get all files from folder               
        files = rdir(folder_list{m,2});
        tier = folder_list{m,1};

        % iterate over all files
        for k = 1:length(files)
            % get filename without extension
            [filepath,name,ext] = fileparts(files(k).name);

            % check if file is in excluded folder
            skip = false;
            for n = 1:length(exclude_folders)
                if startsWith(filepath, exclude_folders{n})
                    skip = true;
                end
            end

            if ~skip && ~isempty(meta.class.fromName(name))
                class_list{id, 1} = name;
                class_list{id, 2} = tier;
                id = id+1;
            end
        end
    end


    % Iterate through all classes identified, and create class inheritance
    % diagram for the particular class
    for k = 1:length(class_list)

        fileBaseName = [class_list{k,2} '_' class_list{k,1} '_INHERITANCE'];
        mytitle = ['Inheritance map for ' class_list{k,1}];
        graphic_file = create_uml_svg({class_list{k,1}}, fileBaseName, mytitle);

        absolute2relative_svglinkpath(graphic_file)     

        if open_svg
            m2uml.display_class_diagram( 'GraphicFile',graphic_file );
        end
    end
end


function graphic_file = create_uml_svg(class_list, fileBaseName, mytitle)

    fqnlist = {};
    cid = 1;
    for k = 1:length(class_list)

        disp(['Processing class: ' class_list{k}])

        % set up list of classes to include in diagram
        fqnlist{cid} = class_list{k};    % Include the class itself
        cid = cid+1;

        tmp = superclasses(class_list{k});  % find all parent classes
        for m = 1:length(tmp)
            % check for known matlab classes we inherit from
            % add more if needed
            % an error will be raised if normal matlab classes are included...
            if ismember(tmp{m},'handle')
                continue;
            elseif ismember(tmp{m},'matlab.mixin.Copyable')
                continue;
            else
                % if it is not a matlab class, add it to the list
                fqnlist{cid} = tmp{m};
                cid = cid+1;
            end
        end
    end

    fqnlist = unique(fqnlist);   % make sure only one of each

    %%  Create a minimal puml-script 
    user.Footer.On          = false;
    user.Header.On          = false;
    user.Title.On           = true;
    user.Diagram.Monospaced = false;
    user.Diagram.Arguments  = false;
    user.Diagram.LinetypeDefault = true;
    user.Property.private   = false;
    user.Property.protected = false;
    user.Method.private     = false;
    user.Method.protected   = false;
    user.TodoFixme.fixme    = false;                     
    user.TodoFixme.todo     = false;                      
    user.TodoFixme.note     = false;

    user.Hyperlink.Class    = true;
    user.Hyperlink.Method   = true;
    user.Hyperlink.Package  = true;
    user.Hyperlink.Property = true;
    user.Hyperlink.Enumeration  = false;
    user.Hyperlink.Event        = false;
    user.Hyperlink.Function     = false;
    user.Hyperlink.TodoFixme    = false;
    %
    user.Tooltip.Class      = true;
    user.Tooltip.Method     = true;
    user.Tooltip.Package    = true;
    user.Tooltip.Property   = true;
    user.Tooltip.Enumeration    = false;
    user.Tooltip.Event          = false;
    user.Tooltip.Function       = false;
    user.Tooltip.TodoFixme      = false;
    %
    user.General.WorkingFolder       = pwd;          % folder of puml and image files
    user.General.PlantUmlExtension   = 'puml';
    user.General.GraphicFormat       = 'svg';        % 'svg', 'png', ...
    user.General.Viewer              = 'web';        % 'web', 'browser' or 'webwindow'
    user.General.PlantUmlJar         = fullfile(pwd, 'plantuml.jar');
    user.General.FileBaseName = fileBaseName;

    puml_script = m2uml.create_PlantUML_script(                     ...
                    'Classes'       , fqnlist                       ...
                ,   'UserOptions'   , user                          ...
                ,   'Title'         , mytitle                       ...
                ,   'FileBaseName'  , fileBaseName ); 

    graphic_file = m2uml.puml2graphic( 'PlantUmlScript',puml_script, 'GraphicFormat','svg' );
end

