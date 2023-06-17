%========================================================================
% CryoGrid MASK class MASK_kml
% selects a region of interest as the inside of a kml track
%
% S. Westermann, Jan 2021
% S. Westermann, Dec 2022
%========================================================================

classdef MASK_kml < matlab.mixin.Copyable

    
    properties
        PARENT
        PARA
        CONST
        STATVAR
    end
    
    methods
        
        function mask = provide_PARA(mask)
            mask.PARA.kml_file_path = [];
            mask.PARA.kml_filename = [];
            mask.PARA.additive = [];
        end

        function mask = provide_STATVAR(mask)

        end
        
        function mask = provide_CONST(mask)
            
        end
        
        function mask = finalize_init(mask)
            
        end
        

        function mask = apply_mask(mask)
            
            [clip_lon, clip_lat] = read_kml(mask, [mask.PARA.kml_file_path mask.PARA.kml_filename]);
            mask_temp = inpolygon(mask.PARENT.STATVAR.latitude,mask.PARENT.STATVAR.longitude, clip_lat,clip_lon);

            if mask.PARA.additive
                mask.PARENT.STATVAR.mask = mask.PARENT.STATVAR.mask | mask_temp;
            else
                mask.PARENT.STATVAR.mask = mask.PARENT.STATVAR.mask & mask_temp;
            end
            
        end
        
        function [Lon, Lat] = read_kml(mask, kmlFile)
            
            [FID msg] = fopen(kmlFile,'rt');
            if FID<0
                error(msg)
            end
            txt = fread(FID,'uint8=>char')';
            fclose(FID);
            expr = '<Placemark.+?>.+?</Placemark>';
            objectStrings = regexp(txt,expr,'match');
            Nos = length(objectStrings);
            for ii = 1:Nos
                % Find Object Name Field
                bucket = regexp(objectStrings{ii},'<name.*?>.+?</name>','match');
                if isempty(bucket)
                    name = 'undefined';
                else
                    % Clip off flags
                    name = regexprep(bucket{1},'<name.*?>\s*','');
                    name = regexprep(name,'\s*</name>','');
                end
                
                % Find Object Description Field
                bucket = regexp(objectStrings{ii},'<description.*?>.+?</description>','match');
                if isempty(bucket)
                    desc = '';
                else
                    % Clip off flags
                    desc = regexprep(bucket{1},'<description.*?>\s*','');
                    desc = regexprep(desc,'\s*</description>','');
                end
                
                geom = 0;
                % Identify Object Type
                if ~isempty(regexp(objectStrings{ii},'<Point', 'once'))
                    geom = 1;
                elseif ~isempty(regexp(objectStrings{ii},'<LineString', 'once'))
                    geom = 2;
                elseif ~isempty(regexp(objectStrings{ii},'<Polygon', 'once'))
                    geom = 3;
                end
                
                switch geom
                    case 1
                        geometry = 'Point';
                    case 2
                        geometry = 'Line';
                    case 3
                        geometry = 'Polygon';
                    otherwise
                        geometry = '';
                end
                
                % Find Coordinate Field
                bucket = regexp(objectStrings{ii},'<coordinates.*?>.+?</coordinates>','match');
                % Clip off flags
                coordStr = regexprep(bucket{1},'<coordinates.*?>(\s+)*','');
                coordStr = regexprep(coordStr,'(\s+)*</coordinates>','');
                % Split coordinate string by commas or white spaces, and convert string
                % to doubles
                coordMat = str2double(regexp(coordStr,'[,\s]+','split'));
                % Rearrange coordinates to form an x-by-3 matrix
                [m,n] = size(coordMat);
                coordMat = reshape(coordMat,3,m*n/3)';
                
                % define polygon in clockwise direction, and terminate
                [Lat, Lon] = poly2ccw(coordMat(:,2),coordMat(:,1));
                if geom==3
                    Lon = [Lon;NaN];
                    Lat = [Lat;NaN];
                end
                
                % Create structure
                kmlStruct(ii).Geometry = geometry;
                kmlStruct(ii).Name = name;
                kmlStruct(ii).Description = desc;
                kmlStruct(ii).Lon = Lon;
                kmlStruct(ii).Lat = Lat;
                kmlStruct(ii).BoundingBox = [[min(Lon) min(Lat);max(Lon) max(Lat)]];
            end
        end
        
%         function [x,y,z] = read_kml(mask, fileName)
%             % READ_KML Reads in (x,y,z) from a GoogleEarth kml file.
%             %
%             %  I have tried to make this code as robust as possible, but it may crash
%             %  or give unexpected resutls if the file is not formatted exactly as
%             %  expected.
%             %
%             % Example:
%             %   [x,y,z] = read_kml('test.kml');
%             %
%             % where test.kml looks like:
%             % <?xml version="1.0" encoding="UTF-8"?>
%             % <kml xmlns="http://earth.google.com/kml/2.1">
%             % <Placemark>
%             % 	<name>test_length</name>
%             % 	<description>junk</description>
%             % 	<LineString>
%             % 		<tessellate>1</tessellate>
%             % 		<coordinates>
%             % -73.65138440596144,40.45517368645169,0 -73.39056199144957,40.52146569128411,0 -73.05890757388369,40.59561213913959,0 -72.80519929505505,40.66961872411046,0 -72.61180114704385,40.72997510603909,0 -72.43718187249095,40.77509309196679,0 </coordinates>
%             % 	</LineString>
%             % </Placemark>
%             % </kml>
%             %
%             % afarris@usgs.gov 2016March09, now can read mulitple sets of coordinates
%             % afarris@usgs.gov 2006November
%             %% open the data file and find the beginning of the data
%             fid=fopen(fileName);
%             if fid < 0
%                 error('could not find file')
%             end
%             % This loop reads the data file one line at a time. If if finds the word
%             % <coordinates>, it knows there is data until it reads the word
%             % </coordinates>.  After loading this data, it keeps reading the file,
%             % looking for another instance of <coordinates> until it finds the word
%             % </kml> which signals that the end of the file has been reached.
%             % Some files have all the data on one line, others have newline charecters
%             % in various points in the file.  I hope this code that works in all cases.
%             done=0;
%             endoffile = 0;
%             ar = 1;
%             while endoffile == 0
%                 while done == 0
%                     junk = fgetl(fid);
%                     f = strfind(junk,'<coordinates>');
%                     ff = strfind(junk,'</kml>');
%                     if ~isempty(f)
%                         done = 1;
%                     elseif  ~isempty(ff)
%                         endoffile = 1;
%                         done = 1;
%                     end
%                 end
%                 if endoffile
%                     break
%                 end
%                 % 'junk' either ends with the word '<coordinates>' OR
%                 % some data follows the word '<coordinates>'
%                 if (f + 13) >= length(junk)
%                     % no data on this line
%                     % done2 is set to zero so the next loop will read the data
%                     done2 = 0;
%                 else
%                     % there is some data in this line following '<coordinates>'
%                     clear f2
%                     f2 = strfind(junk,'</coordinates>');
%                     if ~isempty(f2)
%                         %all data is on this line
%                         % there may be multiple sets of data on this one line
%                         % I read them all
%                         for i = 1 : size(f2,2)
%                             alldata{ar} = junk(f(i)+13:f2(i)-1);
%                             % I add in whitespace b/c sometimes it is missing
%                             alldata{ar+1} = ' ';
%                             ar = ar+2;
%                         end
%                         % done2 is set to one because the next loop does not need to run
%                         done2 = 1;
%                     else
%                         % only some data is on this line
%                         alldata{ar} = junk(f+13:end);
%                         % I add in whitespace b/c sometimes it is missing
%                         alldata{ar+1} = ' ';
%                         ar = ar+2;
%                         % done2 is set to zero so the next loop will read the rest of the data
%                         done2 = 0;
%                     end
%                     % check to see if at end of the file
%                     ff = strfind(junk,'</kml>');
%                     if  ~isempty(ff)
%                         % no more data
%                         endoffile = 1;
%                         break
%                     else
%                         % need to keep looking for more data
%                         done = 0;
%                     end
%                 end
%                 % If not all the data was on the line with the word <coordiate>,
%                 % read in the data
%                 while done2 == 0
%                     % read in line from data file
%                     junk = fgetl(fid);
%                     f = strfind(junk,'</coordinates>');
%                     if isempty(f) == 1
%                         % no ending signal, just add this data to the rest
%                         alldata{ar} = junk;
%                         ar = ar + 1;
%                     else
%                         % ending signal is present
%                         done = 0;
%                         if f < 20
%                             % </coordinates> is in the begining of the line, ergo no data
%                             % on this line; just end the loop
%                             done2 = 1;
%                         else
%                             % the ending signal (</coordinates>) is present: remove it,
%                             % add data to the rest and signal the end of the loop
%                             f2 = strfind(junk,'</coordinates>');
%                             alldata{ar} = junk(1:f2-1);
%                             ar = ar + 1;
%                             done2 = 1;
%                             disp('done with line')
%                         end
%                     end
%                     % check to see if at end of the file
%                     ff = strfind(junk,'</kml>');
%                     if  ~isempty(ff)
%                         % no more data
%                         endoffile = 1;
%                         break
%                     else
%                         % need to keep looking for more data
%                         done = 0;
%                     end
%                 end
%             end
%             fclose(fid);
%             %% get the data into neat vectors
%             %  I have to divide the string into X, Y and Z values.
%             %
%             % This is hard b/c there is no comma between points
%             % (just commans between x and y, and between
%             % y and z)  ie;  -70.0000,42.0000,0 -70.1000,40.10000,0 -70.2,....
%             %
%             % I used to do this by finding commas and spaces, now I use
%             % 'strsplit'!  Thank you Matlab!
%             % 'alldata' is one huge cell
%             % turn alldata into regular vector so it is easier to work with
%             data = cell2mat(alldata);
%             % data is one huge string, split it so there is seperate element for each number
%             C = strsplit(data,{',',' '});
%             % sometimes first and/or last element in C is empty, this causes problems
%             len = size(C,2);
%             if isempty(C{1}) && isempty(C{end})
%                 D = C(2:len-1);
%             elseif isempty(C{1}) && ~isempty(C{end})
%                 D = C(2:end);
%             elseif isempty(C{end}) && ~isempty(C{1})
%                 D = C(1:len-1);
%             end
%             % There has GOT to be a better way to split C into 3 variables!
%             a = 1;
%             for i = 1 : 3: length(D)-2
%                 x(a,1) = str2double(D{i});
%                 a=a+1;
%             end
%             a=1;
%             for i = 2 : 3: length(D)-1
%                 y(a,1) = str2double(D{i});
%                 a=a+1;
%             end
%             a=1;
%             for i = 3 : 3: length(D)
%                 z(a,1) = str2double(D{i});
%                 a=a+1;
%             end
%         end
            
        
        
        %-------------param file generation-----
        function mask = param_file_info(mask)
            mask = provide_PARA(mask);
            
            mask.PARA.STATVAR = [];
            mask.PARA.class_category = 'MASK';
            mask.PARA.options = [];
            
            mask.PARA.comment.kml_file_path = {'folder where kml file is located'};
            mask.PARA.comment.kml_filename = {'kml file name'};
            
            mask.PARA.comment.additive = {'1: region inside kml track added to existing selection; 0: region inside kml track subtracted from existing selection'};
            mask.PARA.default_value.additive = {0};
        end
            
    end
end

