function [s_q, mu_ncv_q, mu_nct_q, mu_st_q, xgrid, ygrid, bbox] = ...
  field_interp(filename, res, fig)

  % run imm on data
  % both data contain two machine worth of data
  c = run_imm(filename);

  c_len = length(c);

  xcoords = cell(1, c_len);
  ycoords = cell(1, c_len);

  for m = 1:length(c_len)
    xcoords{m} = c{m}.x - c{m}.x(1);
    ycoords{m} = c{m}.y - c{m}.y(1);
    s{m} = c{m}.speed;
    mus{m} = c{m}.mu;
  end

  min_X = min([xcoords{:}]);
  max_X = max([xcoords{:}]);
  min_Y = min([ycoords{:}]);
  max_Y = max([ycoords{:}]);

  bbox = struct('min_X', min_X, 'max_X', max_X, ...
    'min_Y', min_Y, 'max_Y', max_Y);

  X = [xcoords{:}];
  Y = [ycoords{:}];
  S = [s{:}];
  MU = vertcat(mus{:});

  [xgrid,ygrid] = meshgrid(bbox.min_X:res:bbox.max_X, ...
    bbox.min_Y:res:bbox.max_Y);

  mu_ncv_q = grid_avg_interp(X, Y, xgrid, ygrid, MU(:,1)');
  mu_nct_q = grid_avg_interp(X, Y, xgrid, ygrid, MU(:,2)');
  mu_st_q = grid_avg_interp(X, Y, xgrid, ygrid, MU(:,3)');
  s_q = grid_avg_interp(X, Y, xgrid, ygrid, S');

  if fig
    figure;
    subplot(1,4,1);
    p1 = imagesc(xgrid(1,:), ygrid(1,:), mu_ncv_q);
    set(p1, 'alphadata', ~isnan(mu_ncv_q));
    set(gca, 'ydir', 'normal', 'yticklabel', [], 'xticklabel', []);
    title('Interpolated p_{ncv}');
    xlabel('Easting');
    ylabel('Northing');
    subplot(1,4,2);
    p2 = imagesc(xgrid(1,:), ygrid(1,:), mu_nct_q);
    set(p2, 'alphadata', ~isnan(mu_nct_q));
    set(gca, 'ydir', 'normal', 'yticklabel', [], 'xticklabel', []);
    title('Interpolated p_{nct}');
    xlabel('Easting');
    ylabel('Northing');
    subplot(1,4,3);
    p3 = imagesc(xgrid(1,:), ygrid(1,:), mu_st_q);
    set(p3, 'alphadata', ~isnan(mu_st_q));
    set(gca, 'ydir', 'normal', 'yticklabel', [], 'xticklabel', []);
    title('Interpolated p_{st}');
    xlabel('Easting');
    ylabel('Northing');
    subplot(1,4,4);
    p4 = imagesc(xgrid(1,:), ygrid(1,:), s_q);
    set(p4, 'alphadata', ~isnan(s_q));
    set(gca, 'ydir', 'normal', 'yticklabel', [], 'xticklabel', []);
    title('Interpolated speed');
    xlabel('Easting');
    ylabel('Northing');
  end

end
