clear all;

year = {'2015', '2017'};
field_name = {'rect1', 'rect2', 'circ'};
path = './mat-files/';

fig_opt = {'Position', [10 10 1900 1800]};

l = 1;
fh = figure(fig_opt{:});
hold on
for m = 1:length(field_name)
  for n = 1:length(year)
    load(strcat(path, field_name{m}, '-', year{n}));
    subplot(3,2,l);
    for p = 1:length(gd)
      ids{p} = gd{p}.id;
      scatter(gd{p}.lon, gd{p}.lat, 8, 'filled');
      xlabel('Longtitude');
      ylabel('Latitude');
      title(strcat(field_name{m}, '-', year{n}));
      set(gca, 'yticklabel', [], 'xticklabel', [], 'fontsize', 15);
      plot_google_map('maptype', 'hybrid');
    end
    legend(unique(ids));
    clear ids;
    l = l + 1;
  end
end
