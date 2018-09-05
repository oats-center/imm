close all;
clear all;

year = {'2015', '2017'};
path = '/backup-disk/data/yaguang-data/'
gps_data_name = '/filesLoadedHistory.mat'
fs_name = '/filesLoadedLocationsFieldShapes.mat'

% load in fieldShape from 2015
load(strcat(path, year{1}, fs_name));
fs = fieldShapes{20};

for m = 1:length(year)
  disp('Year is:');
  disp(year(m));

  gps = load(strcat(path, year{m}, gps_data_name));
  combine_gps = gps.files(gps.fileIndicesCombines); % get only the combine data

  for n = 1:length(combine_gps)
    lonlat = [combine_gps(n).lon combine_gps(n).lat];
    if sum(inShape(fs, lonlat)) > 0
      disp('Combine data number is:');
      disp(n);
      disp('Machine type is:')
      disp(combine_gps(n).id);
      disp('')
      %figure
      %scatter(combine_gps(n).lon(inShape(fs, lonlat)), ...
      %  combine_gps(n).lat(inShape(fs, lonlat)), 5, 'b');
      %axis square
      gps_data = struct('lat', lonlat(inShape(fs, lonlat),2), ...
                        'lon', lonlat(inShape(fs, lonlat),1));
      mname = combine_gps(n).id(~isspace(combine_gps(n).id));
      save(strcat(mname, '-', year{m}, '-', num2str(n)), 'gps_data');
    end
  end
end

% from year 2015, we need:
%   combine data 10
%   combine data 27

% from year 2017, we need:
%   combine data 66
