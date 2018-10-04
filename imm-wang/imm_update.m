function [L, imm_mu_k] = imm_update(imm_mu_km1, z_jbar, S_j, c_j)

  % Get the number of models
  m = length(imm_mu_km1);

  % Compute likelihood for each model
  L = nan(1,m);

  for k = 1:m
    L_log = -0.5 * log(det(2*pi*S_j{k})) - ...
      0.5 * z_jbar{k}' * (S_j{k} \ z_jbar{k});
    L(k) = exp(L_log);
  end

  % Compute normalization constant
  c = sum(L .* c_j);

  % Compute updated model probabilities
  imm_mu_k = L .* c_j / c;
%  if any(imm_mu_k == 0)
%    I = (imm_mu_k == 0);
%    mu_diff = 1e-4;
%    imm_mu_k(~I) = imm_mu_k(~I) - mu_diff;
%    imm_mu_k(I) = 1e-4;
%  end

end %EOF
