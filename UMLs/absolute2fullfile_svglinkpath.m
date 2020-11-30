%==========================================================================
% Converts links in svg files from absolute or relative paths to a relative
% path ('../') that is dynamically prepended with the working directory on
% when the link is pressed. This works without modification on any machine
% independently of the cryogrid install path, but requires the working
% directory of matlab to be set to the UMLs directory of cryogrid in order
% for links to work.
%
% T. Ingemann-Nielsen, November 2020
%==========================================================================

function absolute2fullfile_svglinkpath(graphic_file)
    s=importdata(graphic_file);

    expr_line = "matlab:matlab.desktop.editor.openAndGoToLine\(([^),]+),([^),]+)\)";
    expr_func = "matlab:matlab.desktop.editor.openAndGoToFunction\(([^),]+),([^),]+)\)";

    % iterate over all lines in file
    for m = 1:length(s)
        % find all links to a specific line in specific class
        [startID, endID,tokens, matches] = regexp(s{m},expr_line,'start','end','tokens','match');
        
        % iterate over all links found
        for n = 1:length(matches)
            mstr = matches{n};  % get the path of the file linked to

            pstr = tokens{n}{1}(2:end-1);        % get the link path without quotes
            cs = {pstr; pwd};                    % structure containing link path and current path 
            Schar = char(cs(:));                 % convert to character array
            b = diff(Schar, 1, 1) == 0;          % find difference between the two strings
            cpath = pstr(find(b));               % get the common part of the path

            % Now replace the common part of the path with '..' and
            % regenerate full link information
            mstr2 = strrep(mstr,tokens{n}{1}, ['fullfile(pwd, ''' strrep(pstr, cpath, '..\') ''')']);

            s{m} = strrep(s{m}, mstr, mstr2);  % replace full link in line of file
        end

        [startID, endID,tokens, matches] = regexp(s{m},expr_func,'start','end','tokens','match');

        % iterate over all links found
        for n = 1:length(matches)
            mstr = matches{n};  % get the path of the file linked to

            pstr = tokens{n}{1}(2:end-1);        % get the link path without quotes
            cs = {pstr; pwd};                    % structure containing link path and current path 
            Schar = char(cs(:));                 % convert to character array
            b = diff(Schar, 1, 1) == 0;          % find difference between the two strings
            cpath = pstr(find(b));               % get the common part of the path

            % Now replace the common part of the path with '..' and
            % regenerate full link information
            mstr2 = strrep(mstr,tokens{n}{1}, ['fullfile(pwd, ''' strrep(pstr, cpath, '..\') ''')']);

            s{m} = strrep(s{m}, mstr, mstr2);  % replace full link in line of file
        end

    end


    % Now we just need to save to text file

    fid = fopen(graphic_file,'w'); 
    fprintf(fid, '%s\n',s{:}) ;
    fclose(fid) ;
end