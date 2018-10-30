function [x_0j, P_0j, c_j, mu_ij] = imm_reinit( ...
  x_hat_km1, P_hat_km1, imm_mu_km1, p_ij, ind, dims)

  % Get the number of models
  m = size(x_hat_km1,2);

  % Compute normarlization constants
  c_j = zeros(1,m);
  for jj = 1:m
    for ii = 1:m
      c_j(jj) = c_j(jj) + p_ij(ii,jj) * imm_mu_km1(ii);
    end
  end

  % Compute mixing probabilities
  mu_ij = zeros(m,m);
  for jj = 1:m
    for ii = 1:m
      mu_ij(jj,ii) = p_ij(jj,ii) * imm_mu_km1(jj) / c_j(ii);
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
