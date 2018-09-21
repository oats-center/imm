clear all;

dims = 5;
nmodels = 3;

dt = 1;

a_func = {};
ia_func = {};
a_param = {};

a_func{1} = [];
a_func{2} = @f_turn;
a_func{3} = [];
ia_func{1} = [];
ia_func{2} = @f_turn_inv;
ia_func{3} = [];
a_param{1} = [];
a_param{2} = {dt};
a_param{3} = [];

% Space for model parameters
ind = cell(1,nmodels);
F = cell(1,nmodels);
L = cell(1,nmodels);
Qc = cell(1,nmodels);
A = cell(1,nmodels);
Q = cell(1,nmodels);
H = cell(1,nmodels);

%% IMM %%
% Indices for different state space
ind{1} = [1 2 3 4]';
ind{2} = [1 2 3 4 5]';
ind{3} = [1 2]';

% Model 1
F{1} = [0 0 1 0;
        0 0 0 1;
        0 0 0 0;
        0 0 0 0];

L{1} = [0 0;
        0 0;
        1 0;
        0 1];

q1 = 0.0001;
Qc{1} = q1 * eye(2);

[A{1}, Q{1}] = lti_disc(F{1}, L{1}, Qc{1}, dt);

% Model 2
A{2} = @f_turn_dx;

Qc{2} = 0.01;

L{2} = [0 0 0 0 1]';

Q{2} = L{2} * Qc{2} * L{2}' * dt;

% Model 3
A{3} = [1 0;
        0 1];

Q{3} = [0.0001 0;
        0 0.0001];

% Measurement models
H{1} = [1 0 0 0;
        0 1 0 0];

H{2} = [1 0 0 0 0;
        0 1 0 0 0];

H{3} = [1 0;
        0 1];

% Measurement noise
R = cell(1, nmodels);

r1 = 0.04;
r2 = 0.04;
R{1} = diag([r1 r2]);

r1 = 0.04;
r2 = 0.04;
R{2} = diag([r1 r2]);

r1 = 0.04;
r2 = 0.04;
R{3} = diag([r1 r2]);

% Number of data points
n = 600;

% Space for real states and modes
X_r = zeros(dims,n);
mstate = zeros(1,n);

%%%%%%% Creation of trajectory %%%%%%%
mstate(1:100) = 1;
X_r(:,1) = [0 0 0 0.5 0]';

mstate(101:109) = 2;
%X_r(4,1000) = 0.5;
X_r(5,100) = -deg2rad(20);

mstate(110:209) = 1;
X_r(4,109) = 0.5;

mstate(210:218) = 2;
%X_r(4,2009) = 0.5;
X_r(5,209) = deg2rad(20);

mstate(219:318) = 1;
X_r(4,218) = 0.5;

mstate(319:327) = 2;
%X_r(4,3018) = 0.5;
X_r(5,318) = -deg2rad(20);

mstate(328:427) = 1;
X_r(4,327) = 0.5;

mstate(428:436) = 2;
%X_r(4,4027) = 0.5;
X_r(5,427) = deg2rad(20);

mstate(437:450) = 1;
X_r(4,436) = 0.5;

mstate(451:500) = 3;
X_r(3,450) = 0;
X_r(4,450) = 0;

mstate(501:550) = 1;
X_r(4,500) = 0.5;

mstate(551:600) = 3;
X_r(3,550) = 0;
X_r(4,550) = 0;

%pr_state = zeros(3, length(mstate));
%pr_state(1, find(mstate == 1)) = 1;
%pr_state(2, find(mstate == 2)) = 1;
%pr_state(3, find(mstate == 3)) = 1;

imm_nr_mu = cell(1, 100);
imm_mu = cell(1, 100);

% Initial parameters for IMM
% Initial estimates
x_jhat{1} = [0 0 0 0]';
x_jhat{2} = [0 0 0 0 0]';
x_jhat{3} = [0 0]';

% Initial probabilities
mu_ip = [0.05 0.05 0.9];
%mu_ip = [0.5 0.5];
mu_ip1 = mu_ip;
mu_imm = mu_ip;
a_j = -log(mu_ip1);

% Initial state prediction error covariances
P_jhat{1} = diag([0.1 0.1 0.1 0.1]);
P_jhat{2} = diag([0.1 0.1 0.1 0.1 0.1]);
P_jhat{3} = diag([0.1 0.1]);

% Transition matrix
p_ij = [0.70 0.20 0.10;
        0.20 0.70 0.10;
        0.45 0.45 0.10];
%p_ij = [0.7 0.3;
%        0.3 0.7];

% Generate object state values
for l = 2:n
   st = mstate(l);
   if isstr(a_func{st}) | strcmp(class(a_func{st}),'function_handle')
       X_r(ind{st},l) = feval(a_func{st},X_r(ind{st},l-1),a_param{st});
   else
       X_r(ind{st},l) = A{st}*X_r(ind{st},l-1);
   end
end

% Generate the measurements
Z = zeros(2,n);
for l = 1:n
  Z(:,l) = H{mstate(l)}*X_r(ind{mstate(l)},l) + ...
           gauss_rnd(zeros(size(Z,1),1), R{mstate(l)});
end

% Run IMM on simulated data
for l = 1:n
  % Placeholder for estimated states and labels
  z = zeros(2, length(Z));

  z(1,:) = (Z(1,:) - Z(1,1))'; % x coords
  z(2,:) = (Z(2,:) - Z(2,1))'; % y coords

  z_jbar = cell(1, nmodels);
  S_j = cell(1, nmodels);
  h = cell(1, nmodels);

  % IMM-NR EKF
  % Model-conditioned reinitialization (mixing)
  [x_0j, P_0j, A_j, alpha_j] = imm_nr_reinit( ...
    x_jhat, P_jhat, mu_ip1, p_ij, a_j, ind, dims);

  % Model-conditioned filtering using EKF and regular KF
  for k = 1:nmodels
    [x_jbar{k}, P_jbar{k}, z_jbar{k}, S_j{k}, x_jhat{k}, P_jhat{k}] = ...
      imm_nr_filter( ...
      x_0j{k}(ind{k}), P_0j{k}(ind{k}, ind{k}), A{k}, Q{k}, H{k}, R{k}, ...
      Z(:,l), [], h{k}, a_func{k}, a_param{k});
  end

  % Model probability update
  [mu_ip1, a_j] = imm_nr_update(z_jbar, S_j, A_j, alpha_j);

  % Final combination
  [x_hat, P_hat] = imm_nr_combo(mu_ip1, x_jhat, P_jhat, ind, dims);

  % IMM EKF
  [x_p1, P_p1, c_j1] = eimm_predict(x_jhat, P_jhat, mu_imm, p_ij, ind, ...
                                    dims, A, a_func, a_param, Q);
  [x_ip1, P_ip1, mu_imm, m1, P1] = eimm_update(x_p1, P_p1, c_j1, ind, dims, ...
                                               Z(:,l), H, [], R, []);

  % Save the values
  IMM_NR_x(:,l) = x_hat;
  IMM_NR_P(:,:,l) = P_hat;
  IMM_NR_mu(:,l) = mu_ip1';

  IMM_x(:,l) = m1;
  IMM_P(:,:,l) = P1;
  IMM_mu(:,l) = mu_imm';
end

for l = 1:nmodels
  figure;
  hold on;
  plot(IMM_NR_mu(l,:), 'r', 'linewidth', 1.5);
  plot(IMM_mu(l,:), 'g', 'linewidth', 1.5);
end

figure;
hold on
plot(X_r(1,:), X_r(2,:), 'r', 'linewidth', 1.5);
plot(Z(1,:), Z(2,:), 'g', 'linewidth', 1.5);
plot(IMM_NR_x(1,:), IMM_NR_x(2,:), 'b', 'linewidth', 1.5);
plot(IMM_x(1,:), IMM_x(2,:), 'c', 'linewidth', 1.5);
axis equal
