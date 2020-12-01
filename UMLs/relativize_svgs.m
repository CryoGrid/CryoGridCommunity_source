% localize all relative links in svg files to absolute links relative to
% the current install.

files = rdir('*.svg');

for k = 1:length(files)
    disp(['Converting to relative paths, file: ' files(k).name]);
    absolute2relative_svglinkpath(files(k).name);
end