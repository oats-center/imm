function [mu_ip1] = eimm_update(z_jbar, S_j, A_j, alpha_j)

  % Get the number of models
  m = length(z_jbar);

  % Compute beta for each model and the sum for A_j, beta_j and alpha_j
  beta_j = zeros(1,m);
  aba_sum = zeros(1,m);
  for k = 1:m
    beta_j(k) = 0.5*(z_jbar{k})'*inv(S_j{k})*z_jbar{k} + ...
      0.5*log(norm(2*pi*S_j{k}));
    aba_sum(k) = A_j(k) + beta_j(k) + alpha_j(k);
  end

  % Get the index of minmum sum
  [~, m_I] = min(aba_sum);

  % Compute mu_km
  mu_km = 0;
  for k = 1:m
    mu_km = mu_km + exp((A_j(m_I)+beta_j(m_I)+alpha_j(m_I)) - ...
      (A_j(k)+beta_j(k)+alpha_j(k)));
  end
  mu_km = 1 / mu_km;

  % Compute a_km
  a_km = -log(mu_km);

  % Compute a_kj
  a_kj = zeros(1,m);

  for k = 1:m
    a_kj(k) = A_j(k)+beta_j(k)+alpha_j(k)+a_km- ...
      (A_j(m_I)+beta_j(m_I)+alpha_j(m_I));
    mu_ip1(k) = exp(-a_kj(k));
  end

end %EOF