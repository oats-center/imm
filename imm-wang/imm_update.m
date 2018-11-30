function [L, imm_mu_k] = imm_update(imm_mu_km1, z_jbar, S_j, c_j, l)

  % Get the number of models
  m = length(imm_mu_km1);

  % Compute likelihood for each model
  L = nan(1,m);

  for k = 1:m
    L_log = -0.5 * log(det(2*pi*S_j{k})) - ...
      0.5 * z_jbar{k}' * (S_j{k} \ z_jbar{k});
    L(k) = exp(L_log);
  end

  % Set lowerbound on likelihood values
  % Assumption: NCV model always wins
  if (any(abs(L) < 1e-5)) & (l > 1)
    L(1) = 0.00001;
    L(2) = 0.000001;
  end

  % Compute normalization constant
  c = sum(L .* c_j);

  % Compute updated model probabilities
  imm_mu_k = L .* c_j / c;

end %EOF
