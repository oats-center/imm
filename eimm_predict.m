function [X_p, P_p, A_j, alpha_j] = eimm_predict( ...
  X_ip, P_ip, MU_ip, p_ij, ind, dims, F, a, param, Q)
%IMM_PREDICT  Interacting Multiple Model (IMM) Filter prediction step
%
% Syntax:
%   [X_p, P_p, A_j, alpha_j] = EIMM_PREDICT(
%     X_ip, P_ip, MU_ip, p_ij, ind, dims, A, a, param, Q)
%
% In:
%   X_ip  - Cell array containing N^j x 1 mean state estimate vector for
%           each model j after update step of previous time step
%   P_ip  - Cell array containing N^j x N^j state covariance matrix for
%           each model j after update step of previous time step
%   MU_ip - Vector containing the model probabilities at previous time step
%   p_ij  - Model transition matrix
%   ind   - Indices of state components for each model as a cell array
%   dims  - Total number of different state components in the combined system
%   F     - Dynamic model matrices for each linear model and Jacobians of each
%           non-linear model's measurement model function as a cell array
%   a     - Function handles of dynamic model functions for each model
%           as a cell array
%   param - Parameters of a for each model as a cell array
%   Q     - Process noise matrices for each model as a cell array.
%
% Out:
%   X_p  - Predicted state mean for each model as a cell array
%   P_p  - Predicted state covariance for each model as a cell array
%   A_j  - Minimum of the exponent for odel probabilities in exponent form
%   alpha_j - Parameter needed in the numerically robust IMM
%
% Description:
%   IMM-EKF filter prediction step. If some of the models have linear
%   dynamics standard Kalman filter prediction step is used for those.

  % The number of models
  m = size(X_ip,2);

  % Construct empty cell arrays for ekf_update if a is not specified
  if isempty(a)
      a = cell(1,m);
  end

  % Same for a's parameters
  if isempty(param)
      param = cell(1,m);
  end
%% IMM-NR
  % Compute a_{k-1}^{n}
  a_j = zeros(1,m);
  for l = 1:m
    a_j(m) = -log(MU_ip(l));
  end

  % Compute A_j
  A_j = zeros(1,m);
  for l = 1:m
    A_j(m) = min(-log(MU_ip));
  end

  % Compute alpha_j
  alpha_j = zeros(1,m);
  for l = 1:m
    for k = 1:m
      alpha_j(l) = alpha_j(l) + (-log(m*p_ij(k,l)*exp(a_j(k)-A_j(l))));
    end
  end
%% IMM-NR

  % Normalizing factors for mixing probabilities
%  c_j = zeros(1,m);
%  for l = 1:m
%      for k = 1:m
%          c_j(l) = c_j(l) + p_ij(k,l).*MU_ip(k);
%      end
%  end

  % Mixing probabilities
%  MU_ij = zeros(m,m);
%  for l = 1:m
%      for k = 1:m
%          MU_ij(k,l) = p_ij(k,l) * MU_ip(l) / c_j(k);
%      end
%  end

%% IMM-NR
  % Mixing probabilities
  MU_ij = zeros(m,m);
  for l = 1:m
      for k = 1:m
          MU_ij(k,l) = p_ij(k,l) * exp(-(a_j(k) - A_j(l) + alpha_j(l)));
      end
  end
%% IMM-NR

  % Calculate the mixed state mean for each filter
  X_0j = cell(1,m);
  for l = 1:m
      X_0j{l} = zeros(dims,1);
      for k = 1:m
          X_0j{l}(ind{k}) = X_0j{l}(ind{k}) + X_ip{k}*MU_ij(k,l);
      end
  end

  % Calculate the mixed state covariance for each filter
  P_0j = cell(1,m);
  for l = 1:m
      P_0j{l} = zeros(dims,dims);
      for k = 1:m
          P_0j{l}(ind{k},ind{k}) = P_0j{l}(ind{k},ind{k}) + ...
            MU_ij(k,l)*(P_ip{k} + ...
            (X_ip{k}-X_0j{l}(ind{k}))*(X_ip{k}-X_0j{l}(ind{k}))');
      end
  end

  % Space for predictions
  X_p = cell(1,m);
  P_p = cell(1,m);

  % Make predictions for each model
  for k = 1:m
      [X_p{k}, P_p{k}] = ekf_predict1( ...
        X_0j{k}(ind{k}), P_0j{k}(ind{k},ind{k}), F{k}, Q{k}, a{k}, [], param{k});
  end

end %EOF
