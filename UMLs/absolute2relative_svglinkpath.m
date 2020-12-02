%==========================================================================
% Converts links in svg files from absolute to relative ('../modules') paths
% T. Ingemann-Nielsen, November 2020
%==========================================================================

function absolute2relative_svglinkpath(graphic_file)
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
            mstr2 = strrep(mstr,tokens{n}{1}, ['''' strrep(pstr, cpath, '..\') '''']);

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
            mstr2 = strrep(mstr,tokens{n}{1}, ['''' strrep(pstr, cpath, '..\') '''']);

            s{m} = strrep(s{m}, mstr, mstr2);  % replace full link in line of file
        end

    end


    % Now we just need to save to text file

    fid = fopen(graphic_file,'w'); 
    fprintf(fid, '%s\n',s{:}) ;
    fclose(fid) ;
end