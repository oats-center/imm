function [x_0j, P_0j, A_j, alpha_j] = imm_nr_reinit( ...
  x_ip, P_ip, MU_ip, p_ij, a_j, ind, dims)
%IMM_NR_REINIT Model-conditioned reinitialization step
%
% Syntax:
%   [x_0j, P_0j, A_j, alpha_j] = IMM_NR_REINIT(x_ip, P_ip, MU_ip, p_ij,
%                                              ind, dims)
%
% In:
%   x_ip  - Cell array containing N^j x 1 mean state estimate vector for
%           each model j after update step of previous time step
%   P_ip  - Cell array containing N^j x N^j state covariance matrix for
%           each model j after update step of previous time step
%   MU_ip - Vector containing the model probabilities at previous time step
%   p_ij  - Model transition matrix
%
% Out:
%   x_0j    - Reinitialized (mixed) states
%   P_0j    - Reinitialized (mixed) state covariance
%   A_j     - Minimum of the exponent for model probabilities in exponent form
%   alpha_j - Parameter needed in the numerically robust IMM

  % The number of models
  m = size(x_ip,2);

  % Compute a_{k-1}^{n}
%  a_j = zeros(1,m);
%  for l = 1:m
%    % Hacky trick
%    if MU_ip(l) < 1e-1
%      MU_ip(l) = 1e-1;
%    end
%    a_j(l) = -log(MU_ip(l));
%  end

  % Compute A_j
  A_j = zeros(1,m);
  for l = 1:m
    A_j(l) = min(a_j);
  end

%  %% Verification step (debug)
%  MU_Aj = zeros(1,m);
%  for l = 1:m
%    for k = 1:m
%      MU_Aj(l) = MU_Aj(l) + p_ij(k,l) * MU_ip(l);
%    end
%    MU_Aj(l) = exp(A_j(l))*MU_Aj(l);
%    if MU_Aj(l) <= 0
%      fprintf('\t*************MU_Aj <= 0**************\n');
%    end
%  end
%  alpha_j_temp = -log(MU_Aj);

  % Compute alpha_j
  alpha_j = zeros(1,m);
  for l = 1:m
    for k = 1:m
      alpha_j(l) = alpha_j(l) + p_ij(k,l)*exp(-(a_j(k)-A_j(l)));
    end
    alpha_j(l) = -log(alpha_j(l));
%    if alpha_j(l) < 0
%      fprintf('\t*************alpha_j < 0**************\n');
%    end
  end

  % Mixing probabilities
  MU_ij = zeros(m,m);
  for l = 1:m
    for k = 1:m
        MU_ij(l,k) = p_ij(l,k) * exp(-(a_j(l) - A_j(k) + alpha_j(k)));
    end
  end

  % Calculate the mixed state mean for each filter
  x_0j = cell(1,m);
  for l = 1:m
    x_0j{l} = zeros(dims,1);
    for k = 1:m
        x_0j{l}(ind{k}) = x_0j{l}(ind{k}) + x_ip{k}*MU_ij(k,l);
    end
  end

  % Calculate the mixed state covariance for each filter
  P_0j = cell(1,m);
  for l = 1:m
    P_0j{l} = zeros(dims,dims);
    for k = 1:m
        P_0j{l}(ind{k},ind{k}) = P_0j{l}(ind{k},ind{k}) + ...
          MU_ij(k,l)*(P_ip{k} + ...
          (x_0j{l}(ind{k})-x_ip{k})*(x_0j{l}(ind{k})-x_ip{k})');
    end
  end

end %EOF
