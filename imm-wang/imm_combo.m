function [x_hat, P_hat] = imm_combo(imm_mu_k, x_jhat, P_jhat, ind, dims)

  % Get the number of models
  m = length(x_jhat);

  % Allocate space
  x_hat = zeros(dims,1);
  P_hat = zeros(dims,dims);

  % Output combined state estimate
  for k = 1:m
    x_hat(ind{k}) = x_hat(ind{k}) + x_jhat{k} * imm_mu_k(k);
  end

  % Output combined state covariance
  for k = 1:m
    P_hat(ind{k}, ind{k}) = P_hat(ind{k}, ind{k}) + (P_jhat{k} + ...
      (x_hat(ind{k}) - x_jhat{k}) * (x_hat(ind{k}) - x_jhat{k})') * imm_mu_k(k);
  end
