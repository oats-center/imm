function [x_jbar, P_jbar, z_jbar, S_j, x_jhat, P_jhat, sigma_j] = ...
  imm_nr_lfilter(x_0j, P_0j, A, Q, H, R, Z)

  % Check defaults
  if isempty(A)
  end
  if isempty(Q)
  end
  if isempty(H)
  end

  % Perform filtering
  % Compute x_jbar
  x_jbar = A * x_0j;

  % Compute P_jbar
  P_jbar = A * P_0j * A' + Q;

  % Compute predicted Z
  Z_p = H * x_jbar;

  % Compute innovation
  z_jbar = Z - Z_p;

  % Compute residual covariance
  S_j = R + H * P_jbar * H';

  % Compute sigma_j
  sigma_j = z_jbar' / S_j * z_jbar;
%  if sigma_j > 100
%    P_jbar = A * P_0j * A' + 10 * Q;
%    S_j = R + H * P_jbar * H';
%  end

  % Compute filter gain
  K_j = P_jbar * H' / S_j;

  % Compute updated state
  x_jhat = x_jbar + K_j * z_jbar;

  % Compute updated covariance
  P_jhat = P_jbar - K_j * S_j * K_j';

end %EOF
