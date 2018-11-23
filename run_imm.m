function [gd_mu] = run_imm(m_type, year)
%RUN_IMM perform IMM algorithm on the specified GPS data mat file
%
%   Parameters:
%     m_type - machine type: combine or kart
%     year - data year
%
%  Yang Wang 11/23/2018

  addpath('./imm-wang/');

  path = strcat('./mat-files/', num2str(year));

  if strcmp(m_type, 'combine')
    d_name = '/combine_gd.mat';
  elseif strcmp(m_type, 'kart')
    d_name = '/kart_gd.mat';
  else
    error('Unrecognized machine type!');
  end

  fprintf('The file specified: %s ...\n', strcat(path, d_name));
  load(strcat(path, d_name)); % load in GPS data
  fprintf('Successfully loaded the data!\n')

  % copy the original struct contents
  for m = 1:length(gd)
    for n = 1:length(gd{m})
      for fn = fieldnames(gd{m}{n})'
        gd_mu{m}{n}.(fn{1}) = gd{m}{n}.(fn{1});
      end
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

  % Allocate sensor mean and var
  sensor_mu = zeros(nmodels,1);
  sensor_var = zeros(nmodels,1);

  for mm = 1:length(gd)
    fprintf('\tWorking on dataset %d out of %d total sets:\n', mm, length(gd));
    if length(gd{mm}) == 0
      fprintf('\t\tDataset %d has length 0, move on.\n', mm);
      continue
    end
    for nn = 1:length(gd{mm})
      fprintf('\t\tIMM algorithm is working on %s data ...\n', gd_mu{mm}{nn}.id);

      [x, y] = convert_to_xy(gd_mu{mm}{nn}.lat, gd_mu{mm}{nn}.lon);

      % Add UTM coords to the struct
      gd_mu{mm}{nn}.x = x;
      gd_mu{mm}{nn}.y = y;

      % Placeholder for estimated states and labels
      z = zeros(2, length(gd_mu{mm}{nn}.x));

      z(1,:) = (x - x(1))'; % x coords
      z(2,:) = (y - y(1))'; % y coords

      % IMM-KF-EKF
      imm_mu = mu_i;

      x_jhat{1} = [0 0 0 0]';
      x_jhat{2} = [0 0 0 0 0]';

      P_jhat{1} = diag([0.1 0.1 0.1 0.1]);
      P_jhat{2} = diag([0.1 0.1 0.1 0.1 0.1]);

      outage = 0;
      outage_cnt = 0;

      IMM_LABELS(1) = 1;

      for l = 1:size(z,2)
        % Get the dt from data
        if l > 1
          dt = (gd{mm}{nn}.gpsTime(l) - gd{mm}{nn}.gpsTime(l-1)) / 1000;
          outage = 0;
%          fprintf('\t\tAt iteration %d\n\tx=%f\ty=%f\tdt=%f\n', ...
%            l, z(1,l), z(2,l), dt);
%          fprintf('\t\tm1p=%f\tm2p=%f\n', ...
%            imm_mu(1), imm_mu(2));
          if dt > 5
            fprintf('\t\t**** Possible GPS outage occured! ****\n');
            outage = 1;
            outage_cnt = outage_cnt + 1;
          end
        else
          dt = 1;
          outage = 0;
        end

        % Compute sensor mean and var up to this iteration with labeled
        % sensor data
        [sensor_mu, sensor_var]  = imm_sensor(IMM_LABELS, gd{mm}{nn}.speed(1:l));
%        fprintf('\t\tsensor_mu(1)=%f\tsensor_mu(2)=%f\n', ...
%          sensor_mu(1), sensor_mu(2));
%        fprintf('\t\tsensor_var(1)=%f\tsensor_var(2)=%f\n', ...
%          sensor_var(1), sensor_var(2));
%
        % Update matrices based on sampling interval
        [A{1}, Q{1}] = lti_disc(F{1}, L{1}, Qc{1}, dt);
        Q{2} = L{2} * Qc{2} * L{2}' * dt;
        a_param{2} = {dt};

        % Reinitialize state and covariance (mixing stage)
        [x_0j, P_0j, c_j, ~] = imm_reinit(x_jhat, P_jhat, imm_mu, p_ij, ind, ...
          dims);

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
%       [Li, imm_mu] = imm_update(imm_mu, z_jbar, S_j, c_j);
        [Li, imm_mu] = imm_update(imm_mu, z_jbar, S_j, c_j, ...
          gd{mm}{nn}.speed(l), sensor_mu, sensor_var, outage);

%        fprintf('\t\tLi(1)=%f\tLi(2)=%f\n', ...
%          Li(1), Li(2));

        % Final combination
        [x_hat, P_hat] = imm_combo(imm_mu, x_jhat, P_jhat, ind, dims);

        IMM_X(:,l) = x_hat;
        IMM_P(:,:,l) = P_hat;
        IMM_MU(:,l) = imm_mu';
        [~, IMM_LABELS(l)] = max(imm_mu');
      end

      fprintf('\t\tDone.\n');
      fprintf('\t\tTotal GPS outage count: %d\n', outage_cnt);

      gd_mu{mm}{nn}.mu = IMM_MU';
      gd_mu{mm}{nn}.z = z;
      gd_mu{mm}{nn}.sensor_mu = sensor_mu;
      gd_mu{mm}{nn}.sensor_var = sensor_var;

      clear x;
      clear y;
      clear z;
      clear IMM_X;
      clear IMM_P;
      clear IMM_MU;
      clear IMM_LABELS;

      outage_cnt = 0;
    end
  end

end %EOF
