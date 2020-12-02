% Open index in matlab web browser, to enable links to code
[filepath,name,ext] = fileparts(mfilename('fullpath'));
web('index.html')




function absolute2relative_path(graphic_file)
    s=importdata(graphic_file);

    expr_line = "matlab:matlab.desktop.editor.openAndGoToLine\(([^),]+),([^),]+)\)";
    expr_func = "matlab:matlab.desktop.editor.openAndGoToFunction\(([^),]+),([^),]+)\)";

    for m = 1:length(s)
        [startID, endID,tokens, matches] = regexp(s{m},expr_line,'start','end','tokens','match');
        %disp(['Matches in line: ' num2str(length(matches))]);

        for n = 1:length(matches)
            mstr = matches{n};

            cs = {tokens{n}{1}(2:end-1); pwd};
            Schar = char(cs(:));
            b = diff(Schar, 1, 1) == 0;
            pstr = tokens{n}{1}(2:end-1);
            cpath = pstr(find(b));

            mstr2 = strrep(mstr,tokens{n}{1}, ['fullfile(pwd, ''' strrep(pstr, cpath, '..\') ''')']);

            s{m} = strrep(s{m}, mstr, mstr2);
        end

        [startID, endID,tokens, matches] = regexp(s{m},expr_func,'start','end','tokens','match');
        %disp(['Matches in line: ' num2str(length(matches))]);

        for n = 1:length(matches)
            mstr = matches{n};

            cs = {tokens{n}{1}(2:end-1); pwd};
            Schar = char(cs(:));
            b = diff(Schar, 1, 1) == 0;
            pstr = tokens{n}{1}(2:end-1);
            cpath = pstr(find(b));

            mstr2 = strrep(mstr,tokens{n}{1}, ['fullfile(pwd, ''' strrep(pstr, cpath, '..\') ''')']);

            s{m} = strrep(s{m}, mstr, mstr2);
        end

    end


    % Now we just need to save to text file

    fid = fop
end