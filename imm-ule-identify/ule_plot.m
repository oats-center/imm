function ule_plot(year, save_fig, save_path)

  year = num2str(year);

  fprintf('PLOTTING THE ULE DATA FOR YEAR: %s\n', year);
  if save_fig
    if isempty(save_path)
      error('YOU NEED TO SPECIFY THE OUTPUT PATH IF SAVE_FIG IS SET!\n');
    end
    fprintf('SAVE_FIG = 1\n');
  end

  path = strcat('../mat-files/');
  ule_fname = '/combine_gd_ul.mat';

  load(strcat(path, year, ule_fname)); % load in the ule data, gd_events

  set(0,'defaultTextInterpreter','latex');

  % we don't expect there are more than 4 combines running in the same field
  % simultaneously
  gps_pts_colors = {'r', 'b', 'k', 'c'};

  % unloading events color
  ule_colors = {'g', 'y', 'w', 'm'};

  for m = 1:length(gd_events)
    fprintf('\tIn fs %d:\n', m);
    if length(gd_events{m}) == 0
      fprintf('\t\tNo data in this fs, moving on.\n');
    end
    for n = 1:length(gd_events{m})
      if gd_events{m}{n}.ule_num > 0
        h = figure('pos', [10 10 1500 1500]);
        hold on

        % plot original gps points
        scatter(gd_events{m}{n}.lon, gd_events{m}{n}.lat, 40, ...
          gps_pts_colors{n}, 'filled');

        % check if there are unloading events
        fprintf('\t\t%s has %d unloading events.\n', gd_events{m}{n}.id, ...
          gd_events{m}{n}.ule_num);

        % plot ule in diamonds
        for nn = 1:gd_events{m}{n}.ule_num
          ts = gd_events{m}{n}.ule_ts{nn};
          I = gd_events{m}{n}.gpsTime >= ts(1) & ...
            gd_events{m}{n}.gpsTime <= ts(2);
          scatter(gd_events{m}{n}.lon(I), gd_events{m}{n}.lat(I), 40, ...
            ule_colors{n}, 'd', 'filled');
        end

        % configure the plot
        set(gca, 'xtick', [], 'xticklabel', [], 'ytick', [], ...
          'yticklabel', [], 'fontsize', 30, 'color', 'none');
        id = gd_events{m}{n}.id;
        id(isspace(id)) = [];
        id(regexp(id, '[&]')) = [];
        fig_name = strcat('year\_', year, '\_FS\_', num2str(m), '\_', id, ...
          '\_ULE');
        title(fig_name);
        xlabel('Longitude');
        ylabel('Latitude');

        % configure the legend
        [l, icons] = legend({strcat(id, ' GPS'), strcat(id, ' ULEs')}, ...
          'orientation', 'horizontal', 'location', 'best');
        icons = findobj(icons, '-property', 'Marker', '-and', '-not', ...
          'Marker', 'none');
        set(icons, 'MarkerSize', 10);

        plot_google_map('maptype', 'hybrid');

        pause(2.5)

        if save_fig
          save_name = strcat('year_', year, '_FS_', num2str(m), '_', id, ...
            '_ULE');
          % save and close the figure
          export_fig(strcat(save_path, save_name),  '-painters', ...
            '-transparent');
        end

        close
      else
        fprintf('\t\t%s does not have unloading events.\n', gd_events{m}{n}.id);
      end
    end
  end

end %EOF
