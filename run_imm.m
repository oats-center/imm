function [gd_mu] = run_imm(gd_mat)
%RUN_IMM perform IMM algorithm on the specified GPS data mat file
%
%   Parameters:
%     gd_mat - GPS data file path and name
%
%  Yang Wang 10/25/2018

  addpath('./imm-wang/');
%  addpath('./imm-nr/');

  load(gd_mat); % load in GPS data

  % copy the original struct contents
  for m = 1:length(gd)
    for fn = fieldnames(gd{m})'
      gd_mu{m}.(fn{1}) = gd{m}.(fn{1});
    end
  end

  dims = 5;
  nmodels = 2;

  % Define sampling interval
  dt = nan;

  % Spaces for function handles
  a_func = {};
  a_param = {};

  a_func{1} = [];
  a_func{2} = @f_turn;
  a_param{1} = [];
  a_param{2} = {dt};

  % Space for model parameters
  ind = cell(1,nmodels);
  F = cell(1,nmodels);
  L = cell(1,nmodels);
  Qc = cell(1,nmodels);
  A = cell(1,nmodels);
  Q = cell(1,nmodels);
  R = cell(1,nmodels);
  H = cell(1,nmodels);

  % Indices for indexing different model parameters
  ind{1} = [1 2 3 4]';
  ind{2} = [1 2 3 4 5]';

  % Model 1 - NCV
  F{1} = [0 0 1 0;
          0 0 0 1;
          0 0 0 0;
          0 0 0 0];

  L{1} = [0 0;
          0 0;
          1 0;
          0 1];

  q1 = 1e-2;
  Qc{1} = q1 * eye(2);
%  [A{1}, Q{1}] = lti_disc(F{1}, L{1}, Qc{1}, dt);

  % Model 2 - NCT
  A{2} = @f_turn_dx;

  L{2} = [0 0 0 0 1]';

  Qc{2} = 0.15;

%  Q{2} = L{2} * Qc{2} * L{2}' * dt;

  % Measurement models
  H{1} = [1 0 0 0;
          0 1 0 0];

  H{2} = [1 0 0 0 0;
          0 1 0 0 0];

  % Measurement noise covariance matrices
  R{1} = diag([0.05 0.15]);
  R{2} = diag([0.15 0.45]);

  % Transition matrix
  p_ij = [0.90 0.10;
          0.10 0.90];

  % Initial probabilities
  mu_i = [0.95 0.05];

  for m = 1:length(gd)
    fprintf('\tIMM algorithm is working on %s GPS data\n', gd_mu{m}.id);

    [x, y] = convert_to_xy(gd_mu{m}.lat, gd_mu{m}.lon);

    % Add UTM coords to the struct
    gd_mu{m}.x = x;
    gd_mu{m}.y = y;

    % Placeholder for estimated states and labels
    z = zeros(2, length(gd_mu{m}.x));

    z(1,:) = (x - x(1))'; % x coords
    z(2,:) = (y - y(1))'; % y coords

    % IMM-KF-EKF
    imm_mu = mu_i;

    x_jhat{1} = [0 0 0 0]';
    x_jhat{2} = [0 0 0 0 0]';

    P_jhat{1} = diag([0.1 0.1 0.1 0.1]);
    P_jhat{2} = diag([0.1 0.1 0.1 0.1 0.1]);

    for l = 1:size(z,2)
      % Get the dt from data
      if l < size(z,2)
        dt = (gd{m}.gpsTime(l+1) - gd{m}.gpsTime(l)) / 1000;
%        fprintf('At iteration %d\n\tx=%f\ty=%f\tdt=%f\n', ...
%          l, z(1,l), z(2,l), dt);
%        fprintf('\tm1p=%f\tm2p=%f\n', ...
%          imm_mu(1), imm_mu(2));
%        if dt > 10
%          fprintf('\t**** Possible GPS outage! ****\n');
%        end
      else
        dt = 1;
      end

      % Update matrices based on sampling interval
      [A{1}, Q{1}] = lti_disc(F{1}, L{1}, Qc{1}, dt);
      Q{2} = L{2} * Qc{2} * L{2}' * dt;
      a_param{2} = {dt};

      % Reinitialize state and covariance (mixing stage)
      [x_0j, P_0j, c_j, ~] = imm_reinit(x_jhat, P_jhat, imm_mu, p_ij, ind, dims);

      % Perform model-conditioned filtering
      for k = 1:nmodels
        [x_jbar{k}, P_jbar{k}, z_jbar{k}, ...
          S_j{k}, K_j{k}, x_jhat{k}, P_jhat{k}] = ...
          imm_filter( ...
          x_0j{k}(ind{k}), ...
          P_0j{k}(ind{k}, ind{k}), ...
          A{k}, Q{k}, H{k}, R{k}, z(:,l), a_param{k});
      end

      % Update model probability
      [~, imm_mu] = imm_update(imm_mu, z_jbar, S_j, c_j);

      % Final combination
      [x_hat, P_hat] = imm_combo(imm_mu, x_jhat, P_jhat, ind, dims);

      IMM_X(:,l) = x_hat;
      IMM_P(:,:,l) = P_hat;
      IMM_MU(:,l) = imm_mu';
%{
      % IMM-NR EKF
      % Model-conditioned reinitialization (mixing)
      [x_0j, P_0j, A_j, alpha_j] = imm_nr_reinit( ...
        x_jhat, P_jhat, mu_imm_nr, p_ij, a_j, ind, dims);

      % Model-conditioned filtering
      for k = 1:nmodels
        if k == 3
          [x_jbar{k}, P_jbar{k}, z_jbar{k}, S_j{k}, x_jhat{k}, P_jhat{k}, ...
          sigma_j(k)] = imm_nr_filter( ...
            x_0j{k}(ind{k}), P_0j{k}(ind{k}, ind{k}), A{k}, Q{k}, H{k}, R{k}, ...
            z(:,l), [], [], a_func{k}, a_param{k});
          continue
        end
        [x_jbar{k}, P_jbar{k}, z_jbar{k}, S_j{k}, x_jhat{k}, P_jhat{k}, ~] = ...
          imm_nr_lfilter(x_0j{k}(ind{k}), P_0j{k}(ind{k}, ind{k}), A{k}, ...
          Q{k}, H{k}, R{k}, z(:,l));
      end

      % Model probability update
      [mu_imm_nr, a_j] = imm_nr_update(z_jbar, S_j, A_j, alpha_j);

      % Final combination
      [x_hat, P_hat] = imm_nr_combo(mu_imm_nr, x_jhat, P_jhat, ind, dims);

      % Save the values
      IMM_NR_x(:,l) = x_hat;
      IMM_NR_P(:,:,l) = P_hat;
      IMM_NR_mu(:,l) = mu_imm_nr';
%}
    end

    % add model probabilities to the struct
%    gd_mu{m}.mu = IMM_NR_mu';
%    gd_mu{m}.z = z;
%    gd_mu{m}.xhat = IMM_NR_x;
%    gd_mu{m}.Phat = IMM_NR_P;

    gd_mu{m}.mu = IMM_MU';
    gd_mu{m}.z = z;

%    clear IMM_NR_x;
%    clear IMM_NR_P;
%    clear IMM_NR_mu;
%    clear z_jbar;
%    clear S_j;
    clear x;
    clear y;
    clear z;
    clear IMM_X;
    clear IMM_P;
    clear IMM_MU;

    fprintf('\tDone.\n');
  end

end %EOF
