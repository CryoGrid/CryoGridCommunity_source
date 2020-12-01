%==========================================================================
% Generate index.html for UML files
% T. Ingemann-Nielsen, November 2020
%==========================================================================


filename = 'index.html';

files = rdir('*.svg');

tier1_classes = files(contains({files.name}, 'TIER1'));
tier2_classes = files(contains({files.name}, 'TIER2'));
tier3_classes = files(contains({files.name}, 'TIER3'));
forcing_classes = files(contains({files.name}, 'FORCING'));

%            file_list, Section_name, Description
sections = {{{}, 'TIER 0', 'Tier 0 classes inherit only from intrinsic Matlab classes.'},   % No classes listed, thus empty struct
            {tier1_classes, 'TIER 1',  'Tier 1 classes inherit from Tier 0 classes and are not full functioning classes.'},
            {tier2_classes, 'TIER 2',  ''},
            {tier3_classes, 'TIER 3',  ''},
            {forcing_classes, 'FORCING', ''},
           };

header = {
    '<!DOCTYPE html>'
    '<html>'
    '<head>'
    '<title>CryoGrid Class Ineritance Hierachies</title>'
    '<style>'
    'body {'
    '  background-color: white;'
    '  text-align: left;'
    '  color: black;'
    '  font-family: Arial, Helvetica, sans-serif;'
    '}'
    '</style>'
    '</head>'
    '<body>'
    ''
    '<h1 style="color: #5e9ca0;">CryoGrid Class Ineritance Hierachies</h1>'
    '<p>&nbsp;</p>'};

nbsp = '<p>&nbsp;</p>';
link = '<p><a href="%s" target="_blank" rel="noopener">%s</a></p>';

heading = {
    '<h2 style="color: #2e6c80;"><a id="%s">%s</a></h2>\n'
    '<p>%s</p>\n'};

footer = {
    ''
    '</body>'
    '</html>'};

fid = fopen(filename,'w'); 
fprintf(fid, '%s\n', header{:});

fprintf(fid, '<h2 style="color: #2e6c80;">Contents</h2>\n');
for m = 1:length(sections)
    fprintf(fid, '<p><a href="#%s">Jump to: %s</a></p>\n', sections{m}{2}, sections{m}{2});
end
fprintf(fid, '%s\n', nbsp);

for m = 1:length(sections)
    fprintf(fid, heading{1}, sections{m}{2}, sections{m}{2});
    fprintf(fid, heading{2}, sections{m}{3});
    
    classes = sections{m}{1};
    if ~isempty(classes)
        for k = 1:length(classes)
            [filepath,name,ext] = fileparts(classes(k).name);
            fprintf(fid, link, ['.\' name ext], name);
        end
    end
    fprintf(fid, '%s\n', nbsp);
end

fprintf(fid, '%s\n', footer{:});

fclose(fid) ;