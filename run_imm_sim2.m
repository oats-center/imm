clear all;

addpath('./imm/');
addpath('./imm-nr/');

% Dimensionality of the state space
dims = 4;

% The number of models in use
nmodels = 2;

% Step size
dt = 1;

% Space for function handles and parameters
a_func = {};
ia_func = {};
a_param = {};
h_func = {};
h_param = {};
dh_dx_func = {};

a_func{1} = [];
a_func{2} = [];
ia_func{1} = [];
ia_func{2} = [];
a_param{1} = [];
a_param{2} = [];

% Space for model parameters
ind = cell(1,nmodels);
F   = cell(1,nmodels);
L   = cell(1,nmodels);
Qc  = cell(1,nmodels);
A   = cell(1,nmodels);
Q   = cell(1,nmodels);
H   = cell(1,nmodels);
R   = cell(1,nmodels);

%%% NCV model %%%
ind{1} = [1 2 3 4]';

% Transition matrix for the continous-time velocity model.
F{1} = [0 0 1 0;
        0 0 0 1;
        0 0 0 0;
        0 0 0 0];

% Noise effect matrix for the continous-time system.
L{1} = [0 0;
        0 0;
        1 0;
        0 1];

% Process noise variance
q1 = 1e-4;
Qc{1} = diag([q1 q1]);

% Discretization of the continous-time system.
[A{1}, Q{1}] = lti_disc(F{1},L{1},Qc{1},dt);

%%% ST model %%%
ind{2} = [1 2]';

A{2} = [1 0;
        0 1];

Q{2} = [1e-4 0;
        0 1e-4];

% Measurement models
H{1} = [1 0 0 0;
        0 1 0 0];

H{2} = [1 0;
        0 1];

hdims = 2;

% Noise variances of the measurement models

% NCV model
r1 = 0.01;
r2 = 0.01;
R{1} = diag([r1 r2]);

% ST model
r1 = 0.5;
r2 = 0.5;
R{2} = diag([r1 r2]);

% Number of data points
n = 200;

% Space for real states and modes
X_r = zeros(dims,n);
mstate = zeros(1,n);

%%%%%%% Creation of trajectory %%%%%%%

% Start with constant velocity 1 toward right
mstate(1:40) = 1;
X_r(:,1) = [0 0 1 0]';

% Travel downwards at 40s at velocity 1
mstate(41:90) = 1;
X_r(3,40) = 0;
X_r(4,40) = -1;

% At 90s move stop for 20s
mstate(91:110) = 2;
X_r(3,90) = 0;
X_r(4,90) = 0;

% At 110s commence another turn right with rate -1
mstate(111:160) = 1;
X_r(3,110) = 0;
X_r(4,110) = -1;

% At 16s move straight for 4 seconds
mstate(161:200) = 2;
X_r(3,90) = 0;
X_r(4,90) = 0;

% Generate object state values
for l = 2:n
   st = mstate(l);
   if isstr(a_func{st}) | strcmp(class(a_func{st}),'function_handle')
       X_r(ind{st},l) = feval(a_func{st}, X_r(ind{st},l-1), a_param{st});
   else
       X_r(ind{st},l) = A{st}*X_r(ind{st}, l-1);
   end
end

% Generate the measurements
Z = zeros(hdims,n);
for l = 1:n
  Z(:,l) = H{mstate(l)} * X_r(ind{mstate(l)},l) + ...
    gauss_rnd(zeros(size(Z,1),1), R{mstate(l)});
end

h = plot(X_r(1,:),X_r(2,:),'-g',...
         Z(1,:),Z(2,:),'.',...
         X_r(1,1),X_r(2,1),'ro','linewidth', 1.2, 'MarkerSize', 5);
legend('Real trajectory',...
       'Measurements',...
       'Starting position', 'location', 'southwest');
set(gca,'FontSize', 12);
xlim([-5 95])
ylim([-60 5])

% Initial model probabilities
mu_ip = [0.8 0.2];

% Markov transition probabilities
p_ij = [0.90 0.10;
        0.10 0.90];

%%% Initial estimates %%%

% IMM-EKF
x_ip1{1} = [0 0 1 0]';
x_ip1{2} = [0 0]';
mu_ip1 = mu_ip;

P_ip1{1} = diag([0.1 0.1 0.1 0.1]);
P_ip1{2} = diag([0.1 0.1]);

% Filtering steps
for l = 1:size(Z,2)
  % EKF based IMM
  [x_p1, P_p1, c_j1] = eimm_predict(x_ip1, P_ip1, mu_ip1, p_ij, ind, dims, ...
    A, a_func, a_param, Q);
  [x_ip1, P_ip1, mu_ip1, m1, P1] = imm_update(x_p1, P_p1, c_j1, ind, dims, ...
    Z(:,l), H, R);

  IMM_MM(:,l)   = m1;
  IMM_PP(:,:,l) = P1;
  IMM_MU(:,l)   = mu_ip1';
  IMM_MM_i(:,l) = x_ip1';
  IMM_PP_i(:,l) = P_ip1';
end

% IMM-NR-EKF
x_jhat{1} = [0 0 1 0]';
x_jhat{2} = [0 0]';
a_j = -log(mu_ip);
mu_imm_nr = mu_ip;

P_jhat{1} = diag([0.1 0.1 0.1 0.1]);
P_jhat{2} = diag([0.1 0.1]);

% Filtering steps
for l = 1:size(Z,2)
%  fprintf('\tIteration number: %d\n', l)
  z_jbar = cell(1, nmodels);
  S_j = cell(1, nmodels);

  % IMM-NR EKF
  % Model-conditioned reinitialization (mixing)
  [x_0j, P_0j, A_j, alpha_j] = imm_nr_reinit( ...
    x_jhat, P_jhat, mu_imm_nr, p_ij, a_j, ind, dims);

%  fprintf('P_0j:\n');
%  celldisp(P_0j);
%  fprintf('alpha_j:\n');
%  alpha_j
%  fprintf('A_j:\n');
%  A_j

  % Model-conditioned filtering using EKF and regular KF
  for k = 1:nmodels
    [x_jbar{k}, P_jbar{k}, z_jbar{k}, S_j{k}, x_jhat{k}, P_jhat{k}, ~] = ...
      imm_nr_lfilter(x_0j{k}(ind{k}), P_0j{k}(ind{k}, ind{k}), A{k}, Q{k}, ...
      H{k}, R{k}, Z(:,l));
  end

%  fprintf('x_jbar:\n');
%  celldisp(x_jbar);
%  fprintf('P_jbar:\n');
%  celldisp(P_jbar);
%  fprintf('z_jbar:\n');
%  celldisp(z_jbar);
%  fprintf('x_jhat:\n');
%  celldisp(x_jhat);
%  fprintf('P_jhat:\n');
%  celldisp(P_jhat);
%  fprintf('S_j:\n');
%  celldisp(S_j);

  % Model probability update
  [mu_imm_nr, a_j] = imm_nr_update(z_jbar, S_j, A_j, alpha_j);
%  fprintf('mu_ip1:\n');
%  mu_ip1
%  fprintf('a_j:\n');
%  a_j

  % Final combination
  [x_hat, P_hat] = imm_nr_combo(mu_imm_nr, x_jhat, P_jhat, ind, dims);
%  fprintf('x_hat:\n');
%  x_hat
%  fprintf('P_hat:\n');
%  P_hat
%  fprintf('Z:\n');
%  X_r(:,l)

  % Save the values
  IMM_NR_x(:,l) = x_hat;
  IMM_NR_P(:,:,l) = P_hat;
  IMM_NR_mu(:,l) = mu_imm_nr';
%  IMM_NR_sigma(:,l) = sigma_j';

  % Plot the estimates obtained so far
%  plot(Z(1,1:l),Z(2,1:l),'k.',...
%       IMM_NR_x(1,1:l),IMM_NR_x(2,1:l),'r-',...
%       X_r(1,1:l),X_r(2,1:l),'g-');
%  xlim([-1 7.5])
%  ylim([-3.5 3.5])
%
%  drawnow
end

%%% Plot stuff %%%
figure;
h = plot(X_r(1,:), X_r(2,:), 'g-',...
         Z(1,:), Z(2,:), 'k.', ...
         IMM_MM(1,:), IMM_MM(2,:), '-r', ...
         IMM_NR_x(1,:), IMM_NR_x(2,:), '-b', ...
         'linewidth', 1.2, 'markersize', 10);
legend('True trajectory',...
       'Measurements', ...
       'IMM-EKF', ...
       'IMM-NR-EKF', 'location', 'southwest');
set(gca, 'FontSize', 12);
xlim([-5 95])
ylim([-60 5])

% Determine the real model probabilities
p_models = zeros(nmodels,n);
I1 = find(mstate == 1);
p_models(1,I1) = 1;
I2 = find(mstate == 2);
p_models(2,I2) = 1;

% Plot model 1 probability for each step
figure;
h = plot(1:n, p_models(1,:),'g--',...
         1:n, IMM_MU(1,:)','-r', ...
         1:n, IMM_NR_mu(1,:), '-b', 'linewidth', 1.2);
legend('True',...
      'IMM-EKF', ...
      'IMM-NR-EKF', 'location', 'southwest');
title('Probability of model 1');
ylim([-0.1, 1.1]);
set(gca, 'FontSize', 12);

% Plot model 2 probability for each step
figure;
h = plot(1:n, p_models(2,:),'g--',...
         1:n,IMM_MU(2,:)','-r', ...
         1:n,IMM_NR_mu(2,:), '-b', 'linewidth', 1.2);
legend('True',...
       'IMM-EKF', ...
       'IMM-NR-EKF', 'location', 'northwest');
title('Probability of model 2');
ylim([-0.1, 1.1]);
set(gca, 'FontSize', 12);
