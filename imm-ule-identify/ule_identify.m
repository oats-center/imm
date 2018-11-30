function [combine_gd_ul] = ule_identify(combine_gd_imm, kart_gd_imm)

  % copy the original struct contents
  for m = 1:length(combine_gd_imm)
    for fn = fieldnames(combine_gd_imm{m})'
      combine_gd_ul{m}.(fn{1}) = combine_gd_imm{m}.(fn{1});
    end
  end

  for m = 1:length(combine_gd_ul)
    start_ts = nan;
    end_ts = nan;
    ule_lock = 0;
    ule_num = 0;
    combine_gd_ul{m}.ule_ts = {};

    for mm = 1:length(combine_gd_ul{m}.gpsTime)
      % we assume there is only one cart at any given time ...
      % find the closest matched timestamps
      %TODO: what if there are more than one karts?
      I = find(round(kart_gd_imm{1}.gpsTime / 1000) == ...
        round(combine_gd_ul{m}.gpsTime(mm) / 1000));

      % compute the distance between the combine and the kart
      dis = sqrt((combine_gd_ul{m}.x(mm) - mean(kart_gd_imm{1}.x(I)))^2 + ...
        (combine_gd_ul{m}.y(mm) - mean(kart_gd_imm{1}.y(I)))^2);

      % if the distance is below the hard threshold
      if (dis <= 15) & (abs((combine_gd_ul{m}.mu(mm,1) - ...
        mean(kart_gd_imm{1}.mu(I,1)))) < 0.1)
        if (isnan(start_ts)) & (~ule_lock)
          start_ts = combine_gd_ul{m}.gpsTime(mm);
          ule_lock = 1;
        end
      else
        if (~isnan(start_ts)) & (ule_lock)
          end_ts = combine_gd_ul{m}.gpsTime(mm);
          ule_num = ule_num + 1;
          if (((end_ts - start_ts) / 1000) > 10) | ...
            (((end_ts - start_ts) / 1000) < 300)
            ule_num = ule_num + 1;
            combine_gd_ul{m}.ule_ts{ule_num} = [start_ts end_ts];
          end
          start_ts = nan;
          ule_lock = 0;
        end
      end
    end

    combine_gd_ul{m}.ule_num = ule_num;

    % print info
    fprintf('\tFor %s, total data length %d:\n', ...
      combine_gd_ul{m}.id, length(combine_gd_ul{m}.lat));
    fprintf('\t\tNumber of ULEs: %d\n', combine_gd_ul{m}.ule_num);
  end

  fprintf('\n');

end %EOF
