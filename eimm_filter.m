function [x_jbar, P_jbar, z_jbar, S_j, x_jhat, P_jhat] = eimm_filter( ...
  x_0j, P_0j, A, Q, H, R, z, W, h, V, a, param)

  % Apply defaults
  if isempty(A)
    A = eye(size(x_0j,1));
  end
  if isempty(Q)
    Q = zeros(size(x_0j,1));
  end
  if isempty(W)
    W = eye(size(x_0j,1),size(Q,2));
  end
  if isempty(V)
    V = eye(size(R,1));
  end

  % Perform filtering
  % Evaluate state transition matrix
  if isnumeric(A)
    % nop
  elseif isstr(A) | strcmp(class(A),'function_handle')
    A = feval(A, x_0j, param);
  else
    A = A(x_0j, param);
  end

  % Compute x_jbar
  if isempty(a)
    x_jbar = A * x_0j;
  elseif isnumeric(a)
    x_jbar = a;
  elseif isstr(a) | strcmp(class(a),'function_handle')
    x_jbar = feval(a, x_0j, param);
  else
    x_jbar = a(x_0j, param);
  end

  % Compute P_jbar
  if isnumeric(W)
    % nop
  elseif isstr(W) | strcmp(class(W),'function_handle')
    W = feval(W, x_0j, param);
  else
    W = W(x_0j, param);
  end
  P_jbar = A * P_0j * A' + W * Q * W';
%  P_jbar = A * P_0j * A' + Q;

  % Evaluate measurement matrix
  if isnumeric(H)
    % nop
  elseif isstr(H) | strcmp(class(H),'function_handle')
    H = feval(H, x_jbar, param);
  else
    H = H(x_jbar, param);
  end

  % Compute H*x_jbar
  if isempty(h)
    z_p = H * x_jbar;
  elseif isnumeric(h)
    z_p = h;
  elseif isstr(h) | strcmp(class(h),'function_handle')
    z_p = feval(h, x_jbar, param);
  else
    z_p = h(x_jbar, param);
  end

  if isnumeric(V)
    % nop
  elseif isstr(V) | strcmp(class(V),'function_handle')
    V = feval(V, x_0j, param);
  else
    V = V(x_0j, param);
  end

  % Update everything
  S_j = (V*R*V' + H*P_jbar*H');
%  S_j = R + H*P_jbar*H';
  K_j = P_jbar*H'/S_j;
  z_jbar = z - z_p;
  x_jhat = x_jbar + K_j * z_jbar;
  P_jhat = P_jbar - K_j*S_j*K_j';

end %EOF
