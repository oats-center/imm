function [s_q, mu_ncv_q, mu_nct_q, mu_st_q, xgrid, ygrid, bbox] = ...
  field_interp(filename, res, fig)

  % run imm on data
  % both data contain two machine worth of data
  c = run_imm(filename);

  for m = 1:length(c)
    X = c{m}.x - c{m}.x(1);
    Y = c{m}.y - c{m}.y(1);
    S = c{m}.speed;
    MU = c{m}.mu;

    min_X = min(X);
    max_X = max(X);
    min_Y = min(Y);
    max_Y = max(Y);

    bbox = struct('min_X', min_X, 'max_X', max_X, ...
      'min_Y', min_Y, 'max_Y', max_Y);

    [xgrid,ygrid] = meshgrid(bbox.min_X:res:bbox.max_X, ...
      bbox.min_Y:res:bbox.max_Y);

    mu_ncv_q = grid_avg_interp(X, Y, xgrid, ygrid, MU(:,1)');
    mu_nct_q = grid_avg_interp(X, Y, xgrid, ygrid, MU(:,3)');
    mu_st_q = grid_avg_interp(X, Y, xgrid, ygrid, MU(:,2)');
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

end
