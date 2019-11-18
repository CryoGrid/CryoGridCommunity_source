function zsb = getElevation(Lat, Lon)

load('forcing/subsea6p25kmCellSize_Shelf.mat')    
data = LL6p25km;

zsb = zeros(size(Lat));

for i = 1:length(Lat)
    this_ind = sqrt((data(:,1)-Lat(i)).^2 + (data(:,2)-Lon(i)).^2) == min(sqrt((data(:,1)-Lat(i)).^2 + (data(:,2)-Lon(i)).^2));
    zsb(i) = data(this_ind,3);
end