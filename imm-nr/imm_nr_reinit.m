function [x_0j, P_0j, A_j, alpha_j, mu_ij] = imm_nr_reinit( ...
  x_hat_km1, P_hat_km1, imm_nr_mu_km1, p_ij, a_km1, ind, dims)

  % The number of models
  m = size(x_hat_km1,2);

  % Compute A_j
  A_j = zeros(1,m);
  for jj = 1:m
    A_j(jj) = min(a_km1);
  end

  % Compute alpha_j
  alpha_j = zeros(1,m);
  for jj = 1:m
    for ii = 1:m
      alpha_j(jj) = alpha_j(jj) + p_ij(ii,jj) * exp(-(a_km1(ii)-A_j(jj)));
    end
    alpha_j(jj) = -log(alpha_j(jj));
  end

  % Mixing probabilities
  mu_ij = zeros(m,m);
  for ii = 1:m
    for jj = 1:m
      mu_ij(ii,jj) = p_ij(ii,jj) * exp(-(a_km1(ii) - A_j(jj) + alpha_j(jj)));
    end
  end

  % Reinitialize system state (mixed state)
  x_0j = cell(1,m);
  for k = 1:m
    x_0j{k} = zeros(dims,1);
    for l = 1:m
      x_0j{k}(ind{l}) = x_0j{k}(ind{l}) + x_hat_km1{l} * mu_ij(l,k);
    end
  end

  % Calculate state covariance (mixed covariance)
  P_0j = cell(1,m);
  for k = 1:m
    P_0j{k} = zeros(dims,dims);
    for l = 1:m
      P_0j{k}(ind{l},ind{l}) = P_0j{k}(ind{l},ind{l}) + ...
        mu_ij(l,k) * (P_hat_km1{l} + ...
        (x_hat_km1{l} - x_0j{k}(ind{l})) * (x_hat_km1{l} - x_0j{k}(ind{l}))');
    end
  end

end %EOF
