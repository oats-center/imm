clear all;

set(0,'defaultTextInterpreter','latex');

gd_cb = run_imm('./mat-files/rect1-2015.mat');
gd_tt = run_imm('./mat-files/rect1-2015t.mat');

dlen = length(gd_cb);

I_ul = cell(dlen, 1);

for m = 1:dlen
%  dm = nan(1, length(gd_cb{m}.gpsTime));
  for mm = 1:length(gd_cb{m}.gpsTime)
    I = find(round(gd_tt{1}.gpsTime / 1000) == ...
      round(gd_cb{m}.gpsTime(mm) / 1000));
    dm = sqrt((gd_cb{m}.x(mm) - mean(gd_tt{1}.x(I)))^2 + ...
      (gd_cb{m}.y(mm) - mean(gd_tt{1}.y(I)))^2);
    if dm <= 15
      if (gd_cb{m}.mu(mm,1) > 0.9) & (mean(gd_tt{1}.mu(I,1)) > 0.5)
        I_ul{m} = [I_ul{m} mm];
      end
    end
  end
end

h = figure('pos', [10 10 1500 1500]);
scatter(gd_cb{1}.lon, gd_cb{1}.lat, 40, gd_cb{1}.mu(:,1), ...
  'filled', 'MarkerFaceAlpha', 0.8, 'MarkerEdgeAlpha', 0.8);
c = colorbar;
c.FontSize = 20;
c.Label.Interpreter = 'latex';
caxis([0 1]);
ylabel(c, 'Probability');
set(gca, 'xtick', [], 'xticklabel', [], 'ytick', [], 'yticklabel', [], ...
  'fontsize', 25, 'color', 'none');
plot_google_map('maptype', 'hybrid');
title('Combine 1 NCV Model Probability Map');
xlabel('Longitude');
ylabel('Latitude');
export_fig ./results/rect1_2015_c1pr.png -painters -transparent

h = figure('pos', [10 10 1500 1500]);
scatter(gd_cb{2}.lon, gd_cb{2}.lat, 40, gd_cb{2}.mu(:,1), ...
  'filled', 'MarkerFaceAlpha', 0.8, 'MarkerEdgeAlpha', 0.8);
c = colorbar;
c.FontSize = 20;
c.Label.Interpreter = 'latex';
caxis([0 1]);
ylabel(c, 'Probability');
set(gca, 'xtick', [], 'xticklabel', [], 'ytick', [], 'yticklabel', [], ...
  'fontsize', 25, 'color', 'none');
plot_google_map('maptype', 'hybrid');
title('Combine 2 NCV Model Probability Map');
xlabel('Longitude');
ylabel('Latitude');
export_fig ./results/rect1_2015_c2pr.png -painters -transparent

h = figure('pos', [10 10 1500 1500]);
scatter(gd_tt{1}.lon, gd_tt{1}.lat, 40, gd_tt{1}.mu(:,1), ...
  'filled', 'MarkerFaceAlpha', 0.8, 'MarkerEdgeAlpha', 0.8);
c = colorbar;
c.FontSize = 20;
c.Label.Interpreter = 'latex';
caxis([0 1]);
ylabel(c, 'Probability');
set(gca, 'xtick', [], 'xticklabel', [], 'ytick', [], 'yticklabel', [], ...
  'fontsize', 25, 'color', 'none');
plot_google_map('maptype', 'hybrid');
xlabel('Longitude');
ylabel('Latitude');
title('Tractor NCV Model Probability Map');
export_fig ./results/rect1_2015_tpr.png -painters -transparent

h = figure('pos', [10 10 1500 1500]);
hold on;
scatter(gd_cb{1}.lon, gd_cb{1}.lat, 20, 'r', 'filled');
scatter(gd_tt{1}.lon, gd_tt{1}.lat, 20, 'g', 'filled', ...
  'MarkerFaceAlpha', 0.2, 'MarkerEdgeAlpha', 0.2);
scatter(gd_cb{1}.lon(I_ul{1}), gd_cb{1}.lat(I_ul{1}), 40, 'k^', 'filled');
set(gca, 'xtick', [], 'xticklabel', [], 'ytick', [], 'yticklabel', [], ...
  'fontsize', 25, 'color', 'none');
plot_google_map('maptype', 'hybrid');
legend({'Combine 1 GPS track', 'Grain cart GPS track', ...
  'Combine 1 unloading events'}, 'fontsize', 18, ...
  'orientation', 'horizontal', 'location', 'bestoutside', ...
  'Interpreter','latex');
xlabel('Longitude');
ylabel('Latitude');
title('Combine 1 Unloading Events Map')
export_fig ./results/rect1_2015_ul1.png -painters -transparent

h = figure('pos', [10 10 1500 1500]);
hold on;
scatter(gd_cb{2}.lon, gd_cb{2}.lat, 20, 'b', 'filled', 'MarkerFaceAlpha', 0.2);
scatter(gd_tt{1}.lon, gd_tt{1}.lat, 20, 'g', 'filled', ...
  'MarkerFaceAlpha', 0.2, 'MarkerEdgeAlpha', 0.2);
scatter(gd_cb{2}.lon(I_ul{2}), gd_cb{2}.lat(I_ul{2}), 40, 'k^', 'filled');
set(gca, 'xtick', [], 'xticklabel', [], 'ytick', [], 'yticklabel', [], ...
  'fontsize', 25, 'color', 'none');
plot_google_map('maptype', 'hybrid');
legend({'Combine 2 GPS track', 'Grain cart GPS track', ...
  'Combine 2 unloading events'}, 'fontsize', 18, ...
  'orientation', 'horizontal', 'location', 'bestoutside', ...
  'Interpreter','latex');
xlabel('Longitude');
ylabel('Latitude');
title('Combine 2 Unloading Events Map')
export_fig ./results/rect1_2015_ul2.png -painters -transparent
