clear all;
close all;

map1 = load('pande6088-2015-10-prmap.mat');
map2 = load('pande6130-2015-27-prmap.mat');
map3 = load('pande6130-2017-66-prmap.mat');

figure(1);
scatter(map1.pr_map.coords(:,1), map1.pr_map.coords(:,2), 10, ...
  map1.pr_map.pr(:,1));
axis square;

figure(2);
scatter(map2.pr_map.coords(:,1), map2.pr_map.coords(:,2), 10, ...
  map2.pr_map.pr(:,1));
axis square;

figure(3);
scatter(map3.pr_map.coords(:,1), map3.pr_map.coords(:,2), 10, ...
  map3.pr_map.pr(:,1));
axis square;
