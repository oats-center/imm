function [gd_imm] = run_imm(gd)
%RUN_IMM perform IMM algorithm on the specified GPS data mat file
%
%   Parameters:
%     gd - gps data
%
%  Yang Wang 11/23/2018

  addpath('../imm-wang/');

  % copy the original struct contents
  for m = 1:length(gd)
    for fn = fieldnames(gd{m})'
      gd_imm{m}.(fn{1}) = gd{m}.(fn{1});
    end
  end

  % before running imm, segment the data based on outage
  % find the total outage count and their indices
  for m = 1:length(gd_imm)
    outage_ind = [];
    outage_cnt = 0;
    fprintf('\tFinding outage in dataset %d out of %d total sets:\n', ...
      m, length(gd_imm));

    for l = 1:(length(gd_imm{m}.gpsTime)-1)
      % Get the dt from data
      dt = (gd_imm{m}.gpsTime(l+1) - gd_imm{m}.gpsTime(l)) / 1000;
      if dt > 10
        outage_cnt = outage_cnt + 1;
        outage_ind(outage_cnt) = l;
      end
    end

    % get the outage indices
    gd_imm{m}.outages = outage_ind;
    fprintf('\t\tIn %s, there are %d outages out of %d points.\n', ...
      gd_imm{m}.id, length(gd_imm{m}.outages), length(gd_imm{m}.time));

    if length(gd_imm{m}.outages) == 1
      gd_imm_seg = cell(1, 2);
    elseif length(gd_imm{m}.outages) == 0
      % nop
    else
      gd_imm_seg = cell(1, length(gd_imm{m}.outages));
    end
  end

  fprintf('\n');

  % allocate gps data segments
  for m = 1:length(gd_imm)
    if length(gd_imm{m}.outages) == 1
      fprintf('\tBreak %s data into 2 parts for IMM due to outages:\n', ...
        gd_imm{m}.id);
      cnt = 2;
    elseif length(gd_imm{m}.outages) == 0
      fprintf('\tNo outage for %s data, run normal IMM:\n', ...
       gd_imm{m}.id);
      gd_imm{m} = run_imm_single(gd_imm{m});
      cnt = 0;
    else
      fprintf('\tBreak %s data into %d parts for IMM due to outages:\n', ...
        gd_imm{m}.id, length(gd_imm{m}.outages)+1);
      cnt = length(gd_imm{m}.outages)+1;
    end

    % copy and remove some fields so that structfun can work
    m_type = gd_imm{m}.type;
    id = gd_imm{m}.id;
    outages = gd_imm{m}.outages;
    gd_imm_tmp = rmfield(gd_imm{m}, {'type', 'id', 'outages'});

    if cnt == 2 % when there is only one outage
      % copy everything else
      gd_imm_seg{1} = ...
        structfun(@(x) ...
          x(1:(outages(1))), ...
          gd_imm_tmp, 'uniformoutput', 0);
      gd_imm_seg{1}.id = id;
      gd_imm_seg{1}.type = m_type;

      % run IMM
      fprintf('\t\tRunning segments 1 to %d\n', outages(1));
      gd_imm_seg{1} = run_imm_single(gd_imm_seg{1});

      gd_imm_seg{2} = ...
        structfun(@(x) ...
          x((outages(1)+1):end), ...
          gd_imm_tmp, 'uniformoutput', 0);
      gd_imm_seg{2}.id = id;
      gd_imm_seg{2}.type = m_type;

      % run IMM
      fprintf('\t\tRunning segments %d to end\n', outages(1)+1);
      gd_imm_seg{2} = run_imm_single(gd_imm_seg{2});
    elseif cnt == 0 % when there is no outage
      % nop
    else
      % copy everything else
      gd_imm_seg{1} = ...
        structfun(@(x) ...
          x(1:(outages(1))), ...
          gd_imm_tmp, 'uniformoutput', 0);
      gd_imm_seg{1}.id = id;
      gd_imm_seg{1}.type = m_type;

      % run IMM
      fprintf('\t\tRunning segments 1 to %d\n', outages(1));
      gd_imm_seg{1} = run_imm_single(gd_imm_seg{1});

      for nn = 1:(cnt-2)
        % copy everything else
        gd_imm_seg{nn+1} = ...
          structfun(@(x) ...
            x((outages(nn)+1):(outages(nn+1))), ...
            gd_imm_tmp, 'uniformoutput', 0);
        gd_imm_seg{nn+1}.id = id;
        gd_imm_seg{nn+1}.type = m_type;

        % run IMM
        fprintf('\t\tRunning segments %d to %d\n', ...
          outages(nn)+1, outages(nn+1));
        gd_imm_seg{nn+1} = run_imm_single(gd_imm_seg{nn+1});
      end

      gd_imm_seg{cnt} = ...
        structfun(@(x) ...
          x((outages(cnt-1)+1):end), ...
          gd_imm_tmp, 'uniformoutput', 0);
      gd_imm_seg{cnt}.id = id;
      gd_imm_seg{cnt}.type = m_type;

      % run IMM
      fprintf('\t\tRunning segments %d to end\n', outages(cnt-1)+1);
      gd_imm_seg{cnt} = run_imm_single(gd_imm_seg{cnt});
    end

    if cnt ~= 0
      gd_imm{m}.time = {};
      gd_imm{m}.gpsTime = [];
      gd_imm{m}.lat = [];
      gd_imm{m}.lon = [];
      gd_imm{m}.altitude = [];
      gd_imm{m}.speed = [];
      gd_imm{m}.bearing = [];
      gd_imm{m}.accuracy = [];
      gd_imm{m}.mu = [];
      gd_imm{m}.z = [];
      gd_imm{m}.sensor_mu = [];
      gd_imm{m}.sensor_var = [];
      gd_imm{m}.labels = [];
      gd_imm{m}.x = [];
      gd_imm{m}.y = [];

      % and then we need to put everything back together
      for oo = 1:cnt
        gd_imm{m}.time = [gd_imm{m}.time; gd_imm_seg{oo}.time];
        gd_imm{m}.gpsTime = [gd_imm{m}.gpsTime; gd_imm_seg{oo}.gpsTime];
        gd_imm{m}.lat = [gd_imm{m}.lat; gd_imm_seg{oo}.lat];
        gd_imm{m}.lon = [gd_imm{m}.lon; gd_imm_seg{oo}.lon];
        gd_imm{m}.altitude = [gd_imm{m}.altitude; gd_imm_seg{oo}.altitude];
        gd_imm{m}.speed = [gd_imm{m}.speed; gd_imm_seg{oo}.speed];
        gd_imm{m}.bearing = [gd_imm{m}.bearing; gd_imm_seg{oo}.bearing];
        gd_imm{m}.accuracy = [gd_imm{m}.accuracy; gd_imm_seg{oo}.accuracy];
        gd_imm{m}.mu = [gd_imm{m}.mu; gd_imm_seg{oo}.mu];
        gd_imm{m}.z = [gd_imm{m}.z; gd_imm_seg{oo}.z'];
        gd_imm{m}.sensor_mu = gd_imm_seg{oo}.sensor_mu;
        gd_imm{m}.sensor_var = gd_imm_seg{oo}.sensor_var;
        gd_imm{m}.labels = [gd_imm{m}.labels; gd_imm_seg{oo}.labels];
        gd_imm{m}.x = [gd_imm{m}.x; gd_imm_seg{oo}.x];
        gd_imm{m}.y = [gd_imm{m}.y; gd_imm_seg{oo}.y];
      end
    end
  end
%{
  % IMM parameters
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
    % Setting things up before IMM
    fprintf('\tWorking on dataset %d out of %d total sets:\n', mm, length(gd));
    fprintf('\t\tIMM algorithm is working on %s data ...\n', gd_imm{mm}.id);

    [x, y] = convert_to_xy(gd_imm{mm}.lat, gd_imm{mm}.lon);

    % Add UTM coords to the struct
    gd_imm{mm}.x = x;
    gd_imm{mm}.y = y;

    % Placeholder for estimated states and labels
    z = zeros(2, length(gd_imm{mm}.x));

    z(1,:) = (x - x(1))'; % x coords
    z(2,:) = (y - y(1))'; % y coords

    % IMM-KF-EKF initial conditions
    imm_mu = mu_i;

    x_jhat{1} = [x(1) y(1) 0 0]';
    x_jhat{2} = [x(1) y(1) 0 0 0]';

    P_jhat{1} = diag([0.1 0.1 0.1 0.1]);
    P_jhat{2} = diag([0.1 0.1 0.1 0.1 0.1]);

    outage = 0;
%    outage_cnt = 0;

    IMM_LABELS(1) = 1;

    % IMM starts here
    for l = 1:size(z,2)
      % Get the dt from data
      if l > 1
        dt = (gd_imm{mm}.gpsTime(l) - gd_imm{mm}.gpsTime(l-1)) / 1000;
%        outage = 0;
%        fprintf('\t\tAt iteration %d\n\tx=%f\ty=%f\tdt=%f\n', ...
%          l, z(1,l), z(2,l), dt);
%        fprintf('\t\tm1p=%f\tm2p=%f\n', ...
%          imm_mu(1), imm_mu(2));
%        if dt > 10
%          fprintf('\t\t**** Possible GPS outage occured! ****\n');
%          outage = 1;
%          outage_cnt = outage_cnt + 1;
%        end
      else
        dt = 1;
%        outage = 0;
      end

      % Compute sensor mean and var up to this iteration with labeled
      % sensor data
      [sensor_mu, sensor_var]  = imm_sensor(IMM_LABELS, gd_imm{mm}.speed(1:l));
%      fprintf('\t\tsensor_mu(1)=%f\tsensor_mu(2)=%f\n', ...
%        sensor_mu(1), sensor_mu(2));
%      fprintf('\t\tsensor_var(1)=%f\tsensor_var(2)=%f\n', ...
%        sensor_var(1), sensor_var(2));

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
      [Li, imm_mu] = imm_update(imm_mu, z_jbar, S_j, c_j, ...
        gd_imm{mm}.speed(l), sensor_mu, sensor_var, outage);

%      fprintf('\t\tLi(1)=%f\tLi(2)=%f\n', ...
%        Li(1), Li(2));

      % Final combination
      [x_hat, P_hat] = imm_combo(imm_mu, x_jhat, P_jhat, ind, dims);

      IMM_X(:,l) = x_hat;
      IMM_P(:,:,l) = P_hat;
      IMM_MU(:,l) = imm_mu';
      [~, IMM_LABELS(l)] = max(imm_mu');
    end

    fprintf('\t\tDone.\n');
%    fprintf('\t\tTotal GPS outage count: %d\n', outage_cnt);

    gd_imm{mm}.mu = IMM_MU';
    gd_imm{mm}.z = z;
    gd_imm{mm}.sensor_mu = sensor_mu;
    gd_imm{mm}.sensor_var = sensor_var;
    gd_imm{mm}.labels = IMM_LABELS';

    clear x;
    clear y;
    clear z;
    clear IMM_X;
    clear IMM_P;
    clear IMM_MU;
    clear IMM_LABELS;

%    outage_cnt = 0;
  end
%}
  fprintf('\n');

end %EOF
