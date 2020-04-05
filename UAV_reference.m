%%%%%% Model predictive control 2020 %%%%%
%%% J.M. Schepers & M. van der Leest %%%
%% start script
clear all
clc
% close all

% for j = [1 10 15 20 25 30]     
disp('------------------------------------------------------------------');
disp('          Reference tracking with MPC');
disp('');
disp('------------------------------------------------------------------');
    
%% Load model
load('linearizedmodel.mat')

% Discretize system 
T_pred = 0.1;% sample time
lambda = 0; % time delay
[Ad,Bd,Cd,Dd] = c2dt(sys.A,sys.B,sys.C,T_pred,lambda);
LTI.A= Ad;
LTI.B= Bd;
LTI.C= Cd;
LTI.D= Dd;

OBS_rank = rank(ctrb(Ad',Cd')); %Check observability
rank([eye(size(LTI.A,1))-LTI.A -LTI.B;LTI.C LTI.D]);

%% Initital values
Np = 10; % Prediction horizon
Nc = 5; % Control horiozn
simTime = 6;
T = simTime/T_pred; % sampling steps

% x0= [pi/8 0 0 0 0 0 0.2 0 0.1 0 0 0]'; % <-- state values               % inital conditions x(7) - z states
x0= [pi/8 0 0 0 0 0 0.2 0 0.1 0 0 0]'; % <-- state values X7 outside Xf

PSI_ref = zeros(1,T);
Z_ref = [1*ones(1,T)];
%  Z_ref = [0*ones(1,T/4) 1*ones(1,T/2) 0*ones(1,T/4)];
X_ref = zeros(1,T);
Y_ref = 0*ones(1,T);
r = [PSI_ref;Z_ref;X_ref;Y_ref];   % reference values psi(x5),z(x7),x(x9) and y (x11)

% QR weights
  
Q_weight = 10; % MPC Q state matrix
R_weight = 0.1;  % MPC R control input weight

Q_weight_LQ = 1; % State weight LQ controller
R_weight_LQ = 5; % Input weight LQ controller

disp('------------------------------------------------------------------');
disp('          Simulation time 60 (s)- ');
disp('     x0= [pi/8 0 0 0 0 0 0.2 0 0.1 0 0 0]');
disp('     reference signal Z = 1, X = 0, Y =0');
disp('------------------------------------------------------------------');

%% Initialize state and state estimation vectors
% stating matrices sizes
dim.nx = size(Ad,1); %number of states, 
dim.nu = size(Bd,2); %number of inputs, 
dim.ny = size(Dd,1);
dim.nd = 1;         %number of disturbance inputs
dim.Np = Np;        %prediction horizon
dim.Nc = Nc;        %control horizon
dim.nc = length(Cd(:,1)); % size on c matrix for observability


Q = Q_weight*eye(dim.nx);          % MPC state weight matrix
R = R_weight*eye(dim.nu);            % MPC control input weight matrix


y_ref =r;
x_ref_ots = zeros(dim.nx,T);
u_ref_ots = zeros(dim.nu,T);

G_ref = zeros(dim.nx,dim.nu);
G_ref(5,1) = 1;
G_ref(7,2) = 1;
G_ref(9,3) = 1;
G_ref(11,4) = 1;

x = zeros(length(LTI.A(:,1)),T);    % state trajectory
u = zeros(length(LTI.B(1,:)),T);    % control inputs
y = zeros(length(LTI.C(:,1)),T);    % measurements 
t = zeros(1,T);                 % time vector



Vf = zeros(1,T);                % terminal cost sequence
l = zeros(1,T);                 % stage cost sequence

x(:,1) = x0';

%% LQR
Q_lq = Q_weight_LQ*eye(size(LTI.A));
R_lq = R_weight_LQ*eye(length(Bd(1,:)));
x_lq(:,1) = x0'; % Initial state LQ controller

[G_ref_lq,K_lq,gain_LQ] = LQ_controller(G_ref,LTI,Q_lq,R_lq,T,T_pred);


%% H matrix for MPC
% S= dare(LTI.A,LTI.B,Q,R);  % terminal cost

% adjust terminal cost
S = 1000*eye(dim.nx);

Qbar = blkdiag(kron(eye(dim.Np),Q),S);
Qbar_extnd = Qbar*kron(ones(dim.Np+1,1),eye(dim.nx));
Rbar = kron(R,eye(dim.Np));
Rbar_extnd = Rbar*kron(ones(dim.Np,1),eye(dim.nu));
Sbar = S;

[P,Z] = predmodgen1(LTI,dim);
H = (Z'*Qbar*Z + Rbar);

%% MPC

u_limit_max = 4.9;
u_limit_min = 0.0001;
x_lim_vec = [0.5*pi 1000 0.5*pi 1000 pi 1000 1000 1000 1000 1000 1000 1000]';
x_lim_vec_full = repmat(x_lim_vec,[Np+1 1]);


for k = 1:T
    t(k) = (k-1)*T_pred;
    
    % ~~~~~~~ OTS ~~~~~
    Q_ots = eye(dim.nx);
    R_ots = eye(dim.nu);
    J_ots = blkdiag(Q_ots,R_ots);
    
    A_ots = [eye(dim.nx)-LTI.A -LTI.B; LTI.C zeros(dim.nx,dim.nu)]; %full state information
    b_ots = [zeros(dim.nx,1);G_ref*y_ref(:,k)];
    
    opts = optimoptions('quadprog','Display','off');
    [xr_ur,~,exitflag] = quadprog(J_ots,zeros(dim.nx+dim.nu,1),[],[],A_ots,b_ots,[],[],[],opts);
    
    x_ref_ots(:,k) = xr_ur(1:dim.nx);
    u_ref_ots(:,k) = xr_ur(dim.nx+1:dim.nx+dim.nu);
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
    
    x0 = x(:,k);
%   f = x0'*P'*Qbar*Z - x_ref_ots(:,k)'*P'*Qbar'*Z - u_ref_ots(:,k)'*Rbar_extnd';
    f = x0'*P'*Qbar*Z - x_ref_ots(:,k)'*Qbar_extnd'*Z - u_ref_ots(:,k)'*Rbar_extnd';
    
    % compute control action
    cvx_begin quiet
        variable u_Np(4*Np)
        minimize ( (1/2)*quad_form(u_Np,H) + f*u_Np)
        %input constraints
        u_Np <=  u_limit_max*ones(dim.nu*Np,1);
        u_Np >= -u_limit_max*ones(dim.nu*Np,1);
        % state constraints
        Z(1:Np*dim.nx,:)*u_Np <= -P(1:Np*dim.nx,:)*x0 + x_lim_vec_full(1:Np*dim.nx,:);
        Z(1:Np*dim.nx,:)*u_Np >= -P(1:Np*dim.nx,:)*x0 - x_lim_vec_full(1:Np*dim.nx,:); 
    cvx_end
    
    u(:,k) = u_Np(1:dim.nu); % MPC control action
    
    % apply control action
    x(:,k+1) = Ad*x(:,k) + Bd*u(:,k); 
    y(:,k) = Cd*x(:,k);
    
    %LQ controller
    u_lq(:,k) = -K_lq*x_lq(:,k);
    x_lq(:,k+1) = LTI.A*x_lq(:,k) + LTI.B*u_lq(:,k) + G_ref_lq*r(:,k);
    y_lq(:,k) = LTI.C*x_lq(:,k);
    
    
    % stability analysis    
   [P_s,eigvals,K] = dare(Ad,Bd,Q,R);
    Vf(k) = x(:,k)'*P_s*x(:,k);
    l(k) = 0.5*(x(:,k)'*Q*x(:,k) + u(:,k)'*R*u(:,k));
end

% states_trajectory: Nx12 matrix of trajectory of 12 states
states_trajectory = y';
states_trajectory_lq = y_lq';

% Used to compare different Q and R values
% save_data.Np(j).states_trajectory = states_trajectory;
% save_data.Np(j).states_trajectory_LQ = states_trajectory_lq;
%  end
% %~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Stability analysis
Vf_diff = Vf(2:end) - Vf(1:end-1);
l = l(2:end);
t_s = t(1:end-1);

%%
figure
subplot(2,1,1),stairs(t,states_trajectory(:,7))
hold on
subplot(2,1,1),stairs(t,states_trajectory_lq(:,7));
subplot(2,1,1),stairs(t,y_ref(2,:),'k--')
title('State Trajectory of Z(m)  with MPC','interpreter','latex')
legend('MPC','LQR','Reference trajectory')
hold off

subplot(2,1,2),stairs(t,states_trajectory(:,1))
hold on
subplot(2,1,2),stairs(t,states_trajectory_lq(:,1));
title('state trajectory of $\phi$ [rad]  with MPC','interpreter','latex')
legend('MPC','LQR')


figure
hold on
stairs(t,u(1,:))
stairs(t,u_ref_ots(1,:),'k--')
stairs(t,u_lq(1,:))
title('control input vs u ref')
legend('u MPC','u ref','u LQ')


figure;
hold on;
stairs(t_s,Vf_diff);
stairs(t_s,l);
grid on;
xlabel('time [s]')
% title('state constraints and terminal constraint')
legend('$V_f(f(x,u)) - V_f(x)$','$-l(x,u)$','interpreter','latex');

% %%
% figure
% hold on
% stairs(t,save_data.Np(1).states_trajectory(:,7))
% stairs(t,save_data.Np(10).states_trajectory(:,7))
% stairs(t,save_data.Np(15).states_trajectory(:,7))
% stairs(t,save_data.Np(20).states_trajectory(:,7))
% stairs(t,save_data.Np(25).states_trajectory(:,7))
% stairs(t,save_data.Np(30).states_trajectory(:,7))
% stairs(t,y_ref(2,:),'k--')
% hold off
% legend('Np =1','Np =30','Np =25','Np =20','Np =15','Np =10')
% axis([0 6 0 1.2])
% xlabel('time [s]')
% ylabel('Altitude z[m]')
% grid on
