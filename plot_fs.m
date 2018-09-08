clear all;

year = {'2015', '2016', '2017', '2018'};
path = '/backup-disk/data/yaguang-data/';
fs_name = '/filesLoadedLocationsFieldShapes.mat';

%fs_num_2015 = [20, 57, 21, 58, 1, 2, 6, 42];
%fs_num_2017 = [80, 82, 134, 135, 84, 36, 66, 120, 72, 126];

c = {'r', 'g', 'b', 'y'};

figure;
hold on
for m = 1
  fs = load(strcat(path, year{m}, fs_name));
  for n = 1:length(fs.fieldShapes)
    plot(fs.fieldShapes{n}, 'facecolor', c{m}, 'linestyle', 'none');
    alpha(0.5)
  end
end
axis square
plot_google_map('maptype', 'hybrid')
%{
fs_2015 = load(strcat(path, year{1}, fs_name));
figure;
for m = 1:length(fs_num_2015)
  fprintf('Field shape number is: %d\n', fs_num_2015(m));
  plot(fs_2015.fieldShapes{fs_num_2015(m)}, 'facecolor', 'y', ...
    'linestyle', 'none');
  plot_google_map('maptype', 'hybrid')
  pause
end
%}
