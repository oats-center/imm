function [x_hat, P_hat] = imm_nr_combo(imm_nr_mu, x_jhat, P_jhat, ind, dims)

  % Get the number of models
  m = length(x_jhat);

  % Allocate space
  x_hat = zeros(dims,1);
  P_hat = zeros(dims,dims);

  for k = 1:m
    x_hat(ind{k}) = x_hat(ind{k}) + x_jhat{k} * imm_nr_mu(k);
  end

  for k = 1:m
    P_hat(ind{k}, ind{k}) = P_hat(ind{k}, ind{k}) + (P_jhat{k} + ...
      (x_hat(ind{k}) - x_jhat{k}) * (x_hat(ind{k}) - x_jhat{k})') * ...
      imm_nr_mu(k);
  end
