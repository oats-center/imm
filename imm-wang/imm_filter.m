function [x_jbar, P_jbar, z_jbar, S_j, K_j, x_jhat, P_jhat] = ...
  imm_filter(x_0j, P_0j, A, Q, H, R, Z, param)

  % Check defaults
  if isempty(A)
  end
  if isempty(Q)
  end
  if isempty(H)
  end

  if isstr(A) | strcmp(class(A),'function_handle')
    A = feval(A, x_0j, param);
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
  S_j = H * P_jbar * H' + R;

  % Compute filter gain
  K_j = P_jbar * H' / S_j;

  % Compute updated state
  x_jhat = x_jbar + K_j * z_jbar;

  % Compute updated covariance
  P_jhat = P_jbar - K_j * S_j * K_j';

end %EOF
