function [combine_gd_ul] = ule_identify(combine_gd_imm, kart_gd_imm)

  % allocate unloading events indices
  I_ul = cell(length(combine_gd_imm), 1);

  % copy the original struct contents
  for m = 1:length(combine_gd_imm)
    for fn = fieldnames(combine_gd_imm{m})'
      combine_gd_ul{m}.(fn{1}) = combine_gd_imm{m}.(fn{1});
    end
  end

  for m = 1:length(combine_gd_ul)
    for mm = 1:length(combine_gd_ul{m}.gpsTime)
      % we assume there is only one cart at any given time ...
      % find the closest matched timestamps
      I = find(round(kart_gd_imm{1}.gpsTime / 1000) == ...
        round(combine_gd_ul{m}.gpsTime(mm) / 1000));

      % compute the distance between the combine and the kart
      dis = sqrt((combine_gd_ul{m}.x(mm) - mean(kart_gd_imm{1}.x(I)))^2 + ...
        (combine_gd_ul{m}.y(mm) - mean(kart_gd_imm{1}.y(I)))^2);

      % if the dis is below the hard threshold
      if dis <= 15
        if abs((combine_gd_ul{m}.mu(mm,1) - mean(kart_gd_imm{1}.mu(I,1)))) < 0.1
          I_ul{m} = [I_ul{m} mm];
        end
      end
    end

    % save to struct
    combine_gd_ul{m}.ul = I_ul{m};
    fprintf('FOR %s, TOTAL DATA LENGTH %d:\n', ...
      combine_gd_ul{m}.id, length(combine_gd_ul{m}.lat));
    fprintf('NUMBER OF POINTS IDENTIFIED AS UL EVENTS: %d\n', length(I_ul{m}));
  end

  fprintf('\n');

end %EOF
