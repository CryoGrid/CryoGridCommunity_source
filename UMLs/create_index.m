%==========================================================================
% Generate index.html for UML files
% T. Ingemann-Nielsen, November 2020
%==========================================================================


filename = 'index.html';

files = rdir('*.svg');


tier1_classes = files(contains({files.name}, 'TIER1'));
tier2_classes = files(contains({files.name}, 'TIER2'));
tier3_classes = files(contains({files.name}, 'TIER3'));

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

tier0 = {
    '<h2 style="color: #2e6c80;">TIER 0</h2>'
    '<p>Tier 0 classes inherit only from intrinsic Matlab classes.</p>'};

tier1 = {
    '<h2 style="color: #2e6c80;">TIER 1</h2>'
    '<p>Tier 1 classes inherit from Tier 0 classes and are not full functioning classes.</p>'};

tier2 = {
    '<h2 style="color: #2e6c80;">TIER 2</h2>'};

tier3 = {
    '<h2 style="color: #2e6c80;">TIER 3</h2>'};

footer = {
    ''
    '</body>'
    '</html>'};

fid = fopen(filename,'w'); 
fprintf(fid, '%s\n', header{:});
fprintf(fid, '%s\n', tier0{:});

if ~isempty(tier1_classes)
    for k = 1:length(tier1_classes)
        [filepath,name,ext] = fileparts(tier1_classes(k).name);
        tier1{end+1} = sprintf(link, ['.\' name ext], name);
    end
end
fprintf(fid, '%s\n', tier1{:});
fprintf(fid, '%s\n', nbsp);

if ~isempty(tier2_classes)
    for k = 1:length(tier2_classes)
        [filepath,name,ext] = fileparts(tier2_classes(k).name);
        tier2{end+1} = sprintf(link, ['.\' name ext], name);
    end
end
fprintf(fid, '%s\n', tier2{:});
fprintf(fid, '%s\n', nbsp);

if ~isempty(tier3_classes)
    for k = 1:length(tier3_classes)
        [filepath,name,ext] = fileparts(tier3_classes(k).name);
        tier3{end+1} = sprintf(link, ['.\' name ext], name);
    end
end
fprintf(fid, '%s\n', tier3{:});
fprintf(fid, '%s\n', nbsp);
fclose(fid) ;