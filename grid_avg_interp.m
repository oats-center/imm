function [ avg, xgrid, ygrid, ns ] = grid_avg_interp(x, y, xgrid, ygrid, z)

  res = mode(diff(xgrid(1, :)));

  % Neighborhood radius
  r = res;

  % Chebychev distance -> L∞ norm -> square neighborhoods
  %ns = KDTreeSearcher([x(:), y(:)], 'Distance', 'chebychev');
  ns = KDTreeSearcher([x(:), y(:)], 'Distance', 'euclidean');

  % Do search
  [points, d] = rangesearch(ns, [xgrid(:), ygrid(:)], r);

  % Ignore spot with no points
  IND = find(cellfun(@numel, points) > 0);

  % Calculate weights
  % Modified Shepard's method with squares?
  w = cellfun(@(d) ((r - d) ./ (r .* d)).^2, d(IND), ...
    'UniformOutput', false);

  % Average points in grid
  avg = NaN(size(xgrid));
  avg(IND) = cellfun(@(p, w) w * z(p)' / sum(w), points(IND), w);

end
