% localize all relative links in svg files to absolute links relative to
% the current install.

files = rdir('*.svg');

for k = 1:length(files)
    disp(['Localizing file: ' files(k).name]);
    relative2absolute_svglinkpath(files(k).name);   % Make absolute paths from relative paths (like '../modules/*')
    absolute2absolute_svglinkpath(files(k).name);   % Replace all absolute paths with new absolute paths (to be sure all paths are correct)
end