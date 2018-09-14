function [gd_wpr] = run_imm(gd_mat)
%RUN_IMM perform IMM algorithm on the specified GPS data mat file
%
%   Parameters:
%     gd_mat - GPS data file path and name
%
%  Yang Wang 09/13/2018

  load(gd_mat); % load in GPS data

  % copy the original struct contents
  for m = 1:length(gd)
    for fn = fieldnames(gd{m})'
      gd_wpr{m}.(fn{1}) = gd{m}.(fn{1});
    end
  end

  dims = 5;
  nmodels = 3;

  % Define sampling interval
  dt = NaN;

  % Spaces for function handles
  a_func = {};
  a_param = {};

  a_func{1} = [];
  a_func{2} = [];
  a_func{3} = @f_turn;
  a_param{1} = [];
  a_param{2} = [];
  a_param{3} = {dt};

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
  ind{2} = [1 2]';
  ind{3} = [1 2 3 4 5]';

  % Model 1 - NCV
  F{1} = [0 0 1 0;
          0 0 0 1;
          0 0 0 0;
          0 0 0 0];

  L{1} = [0 0;
          0 0;
          1 0;
          0 1];

  q1 = 0.01;
  Qc{1} = q1 * eye(2);

  % Model 2 - ST
  A{2} = [1 0;
          0 1];

  Q{2} = [0.01 0;
          0 0.01];

  % Model 3 - NCT
  A{3} = @f_turn_dx;

  L{3} = [0 0 0 0 1]';

  Qc{3} = 0.05;

  Q{3} = L{3} * Qc{3} * L{3}' * dt;

  % Measurement models
  H{1} = [1 0 0 0;
          0 1 0 0];

  H{2} = [1 0;
          0 1];

  H{3} = [1 0 0 0 0;
          0 1 0 0 0];

  % Measurement noise covariance matrices
  R{1} = diag([0.1 1]);
  R{2} = diag([1 1]);
  R{3} = diag([0.1 5]);

  % Initial estimates
  x_jhat{1} = [0 0 0 0]';
  x_jhat{2} = [0 0]';
  x_jhat{3} = [0 0 0 0 0]';

  % Initial probabilities
  mu_ip = [0.3 0.4 0.3];
  mu_ip1 = mu_ip;

  % Initial state prediction error covariances
  P_jhat{1} = diag([0.1 0.1 0.1 0.1]);
  P_jhat{2} = diag([0.1 0.1]);
  P_jhat{3} = diag([0.1 0.1 0.1 0.1 0.1]);

  % Transition matrix
  p_ij = [0.70 0.10 0.20;
          0.20 0.10 0.70;
          0.45 0.45 0.10];

  for m = 1:length(gd_wpr)
    fprintf('\tIMM algorithm is working on %s GPS data\n', gd_wpr{m}.id);

    [x, y] = convert_to_xy(gd_wpr{m}.lat, gd_wpr{m}.lon);

    % Add UTM coords to the struct
    gd_wpr{m}.x = x;
    gd_wpr{m}.y = y;

    % Placeholder for estimated states and labels
    z = zeros(2, length(gd_wpr{m}.x));

    z(1,:) = (x - x(1))'; % x coords
    z(2,:) = (y - y(1))'; % y coords

    z_jbar = cell(1, nmodels);
    S_j = cell(1, nmodels);
    h = cell(1, nmodels);

    for l = 1:size(z,2)
      % Get the dt from data
      if l < size(z,2)
        dt = (gd{m}.gpsTime(l+1) - gd{m}.gpsTime(l)) / 1000;
      else
        dt = 1;
      end

%      fprintf('Current iteration:\n');
%      disp(l);

      % Update matrices based on sampling interval
      [A{1}, Q{1}] = lti_disc(F{1}, L{1}, Qc{1}, dt);
      Q{3} = L{3} * Qc{3} * L{3}' * dt;
      a_param{3} = {dt};

      % IMM-NR EKF
      % Model-conditioned reinitialization (mixing)
      [x_0j, P_0j, A_j, alpha_j] = eimm_reinit( ...
        x_jhat, P_jhat, mu_ip1, p_ij, ind, dims);

%      fprintf('Reinitialed conditions:\n');
%      celldisp(x_0j);
%      celldisp(P_0j);

      % Model-conditioned filtering
      for k = 1:nmodels
        [x_jbar{k}, P_jbar{k}, z_jbar{k}, S_j{k}, x_jhat{k}, P_jhat{k}] = ...
          eimm_filter( ...
          x_0j{k}(ind{k}), P_0j{k}(ind{k}, ind{k}), A{k}, Q{k}, H{k}, R{k}, ...
          z(:,l), [], h{k}, [], a_func{k}, a_param{k});
      end

      % Model probability update
      mu_ip1 = eimm_update(z_jbar, S_j, A_j, alpha_j);
      mu_ip1(find(mu_ip1 < 1e-2)) = 1e-2;

%      fprintf('Model probabilities:\n');
%      disp(mu_ip1);

      % Final combination
      [x_hat, P_hat] = eimm_combo(mu_ip1, x_jhat, P_jhat, ind, dims);

%      fprintf('Combined estimates and covariances:\n');
%      disp(x_hat);
%      disp(P_hat);

      EIMM_MM(:,l) = x_hat;
      EIMM_PP(:,:,l) = P_hat;
      EIMM_MU(:,l) = mu_ip1';
    end

    % add model probabilities to the struct
    gd_wpr{m}.mu = EIMM_MU';
    gd_wpr{m}.z = z;
    gd_wpr{m}.xhat = EIMM_MM;

    clear -regexp \w*jbar\b;
    clear -regexp \w*jhat\b;
    clear S_j, z;
    clear h;
    fprintf('\tDone.\n');
  end

end %EOF
