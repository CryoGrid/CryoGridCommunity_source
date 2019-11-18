function Qout = getQ_Davies(lat, lon)
% 
%
%

% updated on May 22,2016 (see "diary May2016.txt"); Paul Overduin

% check that lat and lon are equal in length
if size(lat)~=size(lon)
    disp('lat and lon must be of equal size!')
    return
end
Qout = zeros(size(lat));
 
% load geothermal data
	load forcing/Data_HFip.mat % geothermal heat flux data from Davis et al. (2013, Q3) - data are defined on the same grid as CLIMBER2 SAT and ice sheet data (lon x lat) (240x41)

% find indices
for position = 1:length(lat)
% old version relied on pre-saved indices:
%	load InputData/CLIMBER2/Data_indices_lat_lon_C2  % get index data (lat&lon) for Data matrix
%	[dump, ind_lat(position)] = min(abs(lat_HFip - lat(position)));
%	[dump, ind_lon(position)] = min(abs(lon_HFip - lon(position)));
% new version uses latitude and longitude to interpolate out of HFip
% (QDavies on CLIMBER 2 grid of 240 x 41 cells
% Added reduction to first entry of find, for lat/lon values that are directly in the middle
     indLat = find(abs(lat_HFip-lat(position))==min(abs(lat_HFip-lat(position))), 1, 'first');
     indLon = find(abs(lon_HFip-lon(position))==min(abs(lon_HFip-lon(position))), 1, 'first');
     Qout(position) = HFip(indLon,indLat) / 1000;
end
