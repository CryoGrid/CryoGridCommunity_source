%==========================================================================
% Converts links in svg files from relative ('../') to absolute paths
% according to the current install location.
% T. Ingemann-Nielsen, November 2020
%==========================================================================

function relative2absolute_svglinkpath(graphic_file)
    s=importdata(graphic_file);

    [filepath,name,ext] = fileparts(mfilename('fullpath'));
    pathparts = strsplit(filepath,filesep);
    basepath = strjoin(pathparts(1:end-1), filesep);

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
            
            if startsWith(pstr, '..\')
                % Now replace the .. designator with the full path to the
                % modules directory and regenerate the full link
                % information
                mstr2 = strrep(mstr,tokens{n}{1}, ['''' basepath pstr(3:end) '''']);

                s{m} = strrep(s{m}, mstr, mstr2);  % replace full link in line of file
            end
        end

        [startID, endID,tokens, matches] = regexp(s{m},expr_func,'start','end','tokens','match');

        % iterate over all links found
        for n = 1:length(matches)
            mstr = matches{n};  % get the path of the file linked to

            pstr = tokens{n}{1}(2:end-1);        % get the link path without quotes
            
            if startsWith(pstr, '..\')
                % Now replace the .. designator with the full path to the
                % modules directory and regenerate the full link
                % information
                mstr2 = strrep(mstr,tokens{n}{1}, ['''' basepath pstr(3:end) '''']);

                s{m} = strrep(s{m}, mstr, mstr2);  % replace full link in line of file
            end
        end

    end


    % Now we just need to save to text file

    fid = fopen(graphic_file,'w'); 
    fprintf(fid, '%s\n',s{:}) ;
    fclose(fid) ;
end