function [gd] = get_machine_gd(d_indices, files, fs, fig)
%GET_MACHINE_GD get the GPS data bounded by the specified field shape and
%               machine type
%
%   Parameters:
%     d_indices - indices for different GPS data
%     files - Yaguang's gps data file variable name
%     fs - a particular field shape in enhancedFieldShapes
%     fig - plot toggle
%
%  Yang Wang 11/24/2018

  addpath('../')

  gps = files(d_indices); % get gps by particular indices

  l = 1;
  uniq_gps_data_num = [];
  ids = {};
  % load field shape that fits the field shape
  for n = 1:length(gps)
    if sum(inShape(fs, [gps(n).lon gps(n).lat])) > 0
      fprintf('\tFor this fs, we have the following gps data:\n');
      fprintf('\t\tData number is: %d\n', n);
      fprintf('\t\tMachine type is: %s\n', gps(n).id);
      I = inShape(fs, [gps(n).lon gps(n).lat]);
      %TODO: ASK YAGUANG ABOUT THIS BUG, LAT/LON data length > SPEED, ALT ...
      if length(I) > length(gps(n).altitude)
        I = I(1:length(gps(n).altitude));
      end
      % we only want GPS data that is in the field shape
      gps(n).time = gps(n).time(I);
      gps(n).gpsTime = gps(n).gpsTime(I);
      gps(n).lat = gps(n).lat(I);
      gps(n).lon = gps(n).lon(I);
      gps(n).altitude = gps(n).altitude(I);
      gps(n).speed = gps(n).speed(I);
      gps(n).bearing = gps(n).bearing(I);
      gps(n).accuracy = gps(n).accuracy(I) ;
      % add the unique GPS data num to a placeholder
      if ~(ismember(n, uniq_gps_data_num)) | isempty(uniq_gps_data_num)
        uniq_gps_data_num(:,l) = n;
        ids{l} = gps(n).id;
        l = l + 1;
      end
    end
  end

  fprintf('\n');

  uniq_ids = unique(ids);

  % move on if no gps in this fs
  if length(uniq_ids) == 0
    fprintf('\t\tNo GPS found in this fs!\n');
    fprintf('\n');
    gd = 0;
    return
  end

  fprintf('\t\tThe unique number of machine type is: %d\n', length(uniq_ids));
  fprintf('\t\tThe unique ids: ');
  fprintf('%s | ', uniq_ids{:});
  fprintf('\n');

  if fig
    figure;
    hold on
    % plot the unique gps data
    for mm = 1:length(uniq_gps_data_num)
      scatter(gps(uniq_gps_data_num(mm)).lon, ...
        gps(uniq_gps_data_num(mm)).lat, 8);
    end
    legend
    plot_google_map('maptype', 'hybrid');
  end

  % allocate gps array
  gd = cell(1, length(uniq_ids));

  for mm = 1:length(uniq_ids)
    gd{mm} = struct('type', [], ...
                    'id', [], ...
                    'time', [], ...
                    'gpsTime', [], ...
                    'lat', [], ...
                    'lon', [], ...
                    'altitude', [], ...
                    'speed', [], ...
                    'bearing', [], ...
                    'accuracy', []);
  end

  % we want to concatenate GPS data files with the same id
  for mm = 1:length(uniq_ids)
    for nn = 1:length(uniq_gps_data_num)
      if strcmp(gps(uniq_gps_data_num(nn)).id, uniq_ids{mm})
        gd{mm} = concatenateFiles(gd{mm}, ...
          gps(uniq_gps_data_num(nn)));
      end
    end
  end

  if fig
    figure;
    hold on
    % plot the concatenated unique gps data
    for mm = 1:length(uniq_ids)
      scatter(gd{jj}.lon, gd{mm}.lat, 8);
    end
    legend(uniq_ids);
    plot_google_map('maptype', 'hybrid');
  end

  fprintf('\n');

end %EOF
