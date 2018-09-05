clear all;
close all;

pr_6130_2015 = load('./mat-files/pande6130-2015-27-prmap.mat');
pr_6130_2017 = load('./mat-files/pande6130-2017-66-prmap.mat');
pr_6088_2015 = load('./mat-files/pande6088-2015-10-prmap.mat');

method = {'nearest', 'linear', 'natural'};

x1 = pr_6130_2015.pr_map.coords(:,1);
y1 = pr_6130_2015.pr_map.coords(:,2);
v1 = pr_6130_2015.pr_map.pr(:,1);

x2 = pr_6130_2017.pr_map.coords(:,1);
y2 = pr_6130_2017.pr_map.coords(:,2);
v2 = pr_6130_2017.pr_map.pr(:,1);

x3 = pr_6088_2015.pr_map.coords(:,1);
y3 = pr_6088_2015.pr_map.coords(:,2);
v3 = pr_6088_2015.pr_map.pr(:,1);

[xq1,yq1] = meshgrid(min(x2):1:0, 0:1:max(y2));
vq1 = griddata(x1, y1, v1, xq1, yq1, method{3});

[xq2,yq2] = meshgrid(min(x2):1:0, 0:1:max(y2));
vq2 = griddata(x2, y2, v2, xq2, yq2, method{3});

[xq3,yq3] = meshgrid(min(x2):1:0, 0:1:max(y2));
vq3 = griddata(x3, y3, v3, xq3, yq3, method{3});

figure;
mesh(xq1, yq1, vq1);
view(0,90);
axis square;

figure;
mesh(xq2, yq2, vq2);
view(0,90);
axis square;

figure;
mesh(xq3, yq3, vq3);
view(0,90);
axis square;
