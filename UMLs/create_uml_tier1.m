%==========================================================================
% Generate class inheritance diagrams from main CryoGrid classes Tier 1
% T. Ingemann-Nielsen, November 2020
%==========================================================================

files = rdir('..\modules\TIER_1_processes\**\*.m');
exclude_folders = {'..\modules\TIER_1_processes\INTERACTION'
                   };

% The use of relative paths in matlab is not well supported.
% Carefully test, if adding additional paths above.

process_Tier1 = true;
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


% Now generate one file for all Tier 1 classes
if process_Tier1
    % get all files from folder               
    class_list = {};
    cid = 1;
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
            class_list{cid} = name;
            cid = cid+1;
        end
    end

    fileBaseName = 'TIER1_CLASS_INHERITANCE';
    mytitle = ['Inheritance map for TIER 1 classes'];
    graphic_file = create_uml_svg(class_list, fileBaseName, mytitle);

    absolute2relative_svglinkpath(graphic_file)     

    if open_svg
        m2uml.display_class_diagram( 'GraphicFile',graphic_file );
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

