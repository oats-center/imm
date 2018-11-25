function [gd_events] = imm_ule_identify(year)

  fprintf('IMM_ULE_IDENTIFY algorithm started ...\n');

  % obtain the correct paths and filenames
  year = num2str(year);
  if str2num(year) == 2016
    year = strcat(year, '_SynchedAccordingToAnimiations');
  end
  path = strcat('/backup-disk/GpsTrackData/', year);
  gps_fname = '/filesLoadedHistory.mat';
  fs_fname = '/enhancedFieldShapes.mat';

  fprintf('YEAR: %s\n', year);

  load(strcat(path, gps_fname)); % load in files
  load(strcat(path, fs_fname)); % load in enhancedFieldShapes

  combine_gd = files(fileIndicesCombines); % get the combine data
  kart_gd = files(fileIndicesGrainKarts); % get the grain kart data

  % allocate cell array
  gd_events = cell(1, length(enhancedFieldShapes));

  % go through each field shape
  %   1. get gps data within a field shape
  %   2. perform IMM
  %   3. identify unloading/loading events
  %   4. save to cell array
  for m = 1:length(enhancedFieldShapes)
    fprintf('ON FS %d OUT OF %d TOTAL FS:\n\n', m, length(enhancedFieldShapes));
    % obtain the combine, kart gps data within a fs
    fprintf('GETTING COMBINE DATA ...\n')
    combine_gd = get_machine_gd(fileIndicesCombines, files, ...
      enhancedFieldShapes{m}, 0);

    fprintf('GETTING KART DATA ... \n')
    kart_gd = get_machine_gd(fileIndicesGrainKarts, files, ...
      enhancedFieldShapes{m}, 0);

    % in case one of them does not have data, do not process, go on to the next
    if (~iscell(combine_gd)) | (~iscell(kart_gd))
      fprintf('ONE OR MORE DOES NOT HAVE GPS IN THIS FS, MOVE ON TO NEXT!\n');
      fprintf('\n');
      continue
    end

    % perform IMM on both data sets
    fprintf('RUNNING IMM ON COMBINE DATA ...\n');
    combine_gd_imm = run_imm(combine_gd);

    fprintf('RUNNING IMM ON KART DATA ...\n');
    kart_gd_imm = run_imm(kart_gd);

    % identify the unloading/loading events
    fprintf('IDENTIFYING UNLOAD EVENTS ...\n')
    combine_gd_ul = ule_identify(combine_gd_imm, kart_gd_imm);

    gd_events{m} = combine_gd_ul;
  end

  fprintf('IMM_ULE_IDENTIFY algorithm finished!\n');

end %EOF
