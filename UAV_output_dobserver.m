%%%%%% Model predictive control 2020 %%%%%
%%% J.M. Schepers & M. van der Leest %%%
%% start script
clear all
clc
close all

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

% ~~~~~~~~~~~~ Obsever not all states can be measured 
C = [0 0 0 0 1 0 0 0 0 0 0 0; % state psi 
     0 0 0 0 0 0 1 0 0 0 0 0; % state z  
     0 0 0 0 0 0 0 0 1 0 0 0; % state x
     0 0 0 0 0 0 0 0 0 0 1 0]; % state y
%~~~~~~~~~~~~~~~~~~~~~~~~~

OBS_rank = rank(ctrb(Ad',C')); %Check observability
rank([eye(size(LTI.A,1))-LTI.A -LTI.B;LTI.C LTI.D]);

%% Initital values
Np =10; % Prediction horizon
Nc = 12; % Control horiozn
simTime = 6;
T = simTime/T_pred; % sampling steps

x0= [pi/4 0 0 0 pi/8 0 0.2 0 1 0 0.3 0]';                  % inital conditions x(7) - z states
xtilde0 = zeros(1,12);                              %observer initial state

PSI_ref = zeros(1,T);
% Z_ref = 1*ones(1,T);
Z_ref = [0*ones(1,T/4) 1*ones(1,T/2) 0*ones(1,T/4)];
X_ref = zeros(1,T);
Y_ref = zeros(1,T);
r = [PSI_ref;Z_ref;X_ref;Y_ref];   % reference values psi(x5),z(x7),x(x9) and y (x11)
% r = [PSI_ref;Z_ref];    % Onlye reference signal on psi and Z



% QR weights
Q_weight = 10; % MPC Q state matrix
R_weight = 0.1;  % MPC R control input weight

Q_kf_weight =1; % Observer weight state matrix
R_kf_weight =1.05; % Observer weight control input matrix

Q_weight_LQ = 1; % weight LQ controller states
R_weight_LQ = 5; % weight LQ controller inputs

%% Initialize state and state estimation vectors
% stating matrices sizes
dim.nx = size(Ad,1); %number of states, 
dim.nu = size(Bd,2); %number of inputs, 
dim.ny = size(Dd,1);
dim.nd = 1; %number of disturbance inputs
dim.Np = Np; %prediction horizon
dim.Nc = Nc; %control horizon
dim.nc = length(C(:,1)); % size on c matrix for observability
dim.r = size(r,1);

Q = Q_weight*eye(dim.nx);          % MPC state weight matrix
R = R_weight*eye(dim.nu);            % MPC control input weight matrix

Q_KF = Q_kf_weight*eye(dim.nx);         % Obsever Q weight
R_KF = R_kf_weight*eye(dim.nc);    % Obsever R input weight


y_ref =r;
x_ref_ots = zeros(dim.nx,T);
u_ref_ots = zeros(dim.nu,T);

G_ref = zeros(dim.nx,dim.r);
% gref for only ref 2 states psi and z
% G_ref(5,1) = 1;
% G_ref(7,2) = 1;

% for G_ref 12 by 4
G_ref(5,1) = 1;
G_ref(7,2) = 1;
G_ref(9,3) = 1;
G_ref(11,4) = 1;

xtilde = zeros(dim.nx,T);       %estimated states
ytilde = zeros(dim.nc,T);       % measurements 

x = zeros(length(Ad(:,1)),T);    % state trajectory
u = zeros(length(Bd(1,:)),T);    % control inputs
y = zeros(length(C(:,1)),T);    % measurements 
t = zeros(1,T);                 % time vector

Vf = zeros(1,T);                % terminal cost sequence
l = zeros(1,T);                 % stage cost sequence

x(:,1) = x0';                   % initial state original system
xtilde(:,1) = xtilde0';         % initial state observer
x_lq(:,1) = xtilde0';                % Initial state LQ controller
%% Observer
[~,obs_eig,K] = dare(Ad',C',Q_KF,R_KF);
K = K';


%% H matrix for MPC
S=dare(LTI.A,LTI.B,Q,R); % terminal cost

S = 1000*eye(dim.nx);

Qbar = blkdiag(kron(eye(dim.Np),Q),S);
Qbar_extnd = Qbar*kron(ones(dim.Np+1,1),eye(dim.nx));
Rbar = kron(R,eye(dim.Np));
Rbar_extnd = Rbar*kron(ones(dim.Np,1),eye(dim.nu));
Sbar = S;

[P,Z] = predmodgen1(LTI,dim);
             
H = (Z'*Qbar*Z + Rbar);

%% LQR
Q_lq = Q_weight_LQ*eye(size(LTI.A));
R_lq = R_weight_LQ*eye(length(LTI.B(1,:)));

% two reference states psi and Z
% [G_ref_lq,K_lq] = LQ_controller_ref2(G_ref,LTI,Q_lq,R_lq,T,T_pred);

% Four reference states [Psi, Z, X Y]
[G_ref_lq,K_lq] = LQ_controller(G_ref,LTI,Q_lq,R_lq,T,T_pred);

%% MPC

u_limit_max = 4.9;
u_limit_min = 0.0001;
x_lim_vec = [0.5*pi 1000 0.5*pi 1000 pi 1000 1000 1000 1000 1000 1000 1000]';
x_lim_vec_full = repmat(x_lim_vec,[Np 1]);


for k = 1:T
    t(k) = (k-1)*T_pred;
    
    %%%%%%%% OTS %%%%%%%%%%
    Q_ots = eye(dim.nx);
    R_ots = eye(dim.nu);
    J_ots = blkdiag(Q_ots,R_ots);
    
    A_ots = [eye(dim.nx)-Ad -Bd; Cd zeros(size(Cd,1),dim.nu)];%full state information
    b_ots = [zeros(dim.nx,1);G_ref*y_ref(:,k)];
    
    opts = optimoptions('quadprog','Display','off');
    [xr_ur,~,exitflag] = quadprog(J_ots,zeros(dim.nx+dim.nu,1),[],[],A_ots,b_ots,[],[],[],opts);
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    x_ref_ots(:,k) = xr_ur(1:dim.nx);
    u_ref_ots(:,k) = xr_ur(dim.nx+1:dim.nx+dim.nu);

    % computing f
    x0 = x(1:dim.nx,k); %Estimated state
    f = x0'*P'*Qbar*Z - x_ref_ots(:,k)'*Qbar_extnd'*Z - u_ref_ots(:,k)'*Rbar_extnd';  % Original
%   f = x0'*P'*Qbar*Z - x_ref_ots(:,k)'*P'*Qbar'*Z - u_ref_ots(:,k)'*Rbar_extnd';     % With addditional P value for Xref
    
    % compute control action
    cvx_begin quiet
        variable u_Np(4*Np)
        minimize ( (1/2)*quad_form(u_Np,H) + f*u_Np)
        %input constraints
        u_Np <=  u_limit_max*ones(dim.nu*Np,1);
        u_Np >=  -u_limit_max*ones(dim.nu*Np,1);
        % state constraints
        Z(1:Np*dim.nx,:)*u_Np <= -P(1:Np*dim.nx,:)*x0 + x_lim_vec_full;
        Z(1:Np*dim.nx,:)*u_Np >= -P(1:Np*dim.nx,:)*x0 - x_lim_vec_full; 
    cvx_end
    
    u(:,k) = u_Np(1:dim.nu); % MPC control action
 
    % apply control action
    x(:,k+1) = LTI.A*x(:,k) + LTI.B*u(:,k); %
    y(:,k) = C*x(:,k); % 

   %Augmented observer
   
    ytilde(:,k) = C*xtilde(:,k);
    xtilde(:,k+1) = LTI.A*xtilde(:,k) + LTI.B*u(:,k) + K*(y(:,k) - ytilde(:,k));
   
    %LQ controller
    u_lq(:,k) = -K_lq*x_lq(:,k);
    x_lq(:,k+1) = Ad*x_lq(:,k) + Bd*u_lq(:,k) + G_ref_lq*r(:,k);
    y_lq(:,k) = C*x_lq(:,k);
    
end

% states_trajectory: Nx12 matrix of trajectory of 12 states
states_trajectory_est = ytilde';
states_trajectory = y';
states_trajectory_lq = y_lq';

%%
figure
subplot(2,1,1),stairs(t,states_trajectory_est(:,2))
hold on
subplot(2,1,1),stairs(t,states_trajectory_lq(:,2))
subplot(2,1,1),stairs(t,x_ref_ots(7,:),'k--')
legend('altitude MPC observer z[m]','altitude LQ z[m]','reference altitude z[m]')
hold off

subplot(2,1,2),stairs(t,states_trajectory_est(:,3))
hold on
subplot(2,1,2),stairs(t,states_trajectory_est(:,4))
subplot(2,1,2),stairs(t,states_trajectory_est(:,1))
legend('MPC x[m]','MPC y[m]','MPC $\psi[rad]$','interpreter','latex')
hold off

figure
stairs(t,states_trajectory_est(:,2))
hold on
stairs(t,states_trajectory(:,2))
stairs(t,x_ref_ots(7,:),'k--')
legend('$\tilde{Z} $ [m]', 'Z[m]','interpreter','latex')
axis([0 6 0 1.2])
xlabel('Time [s]')
ylabel('Altitude z[m]')
% title('MPC obsever vs real state')

figure
stairs(t,u(1,:))
hold on
stairs(t,u_ref_ots(1,:),'k--')
title('Control input u1 vs reference')



