%%%%%% Model predictive control 2020 %%%%%
%%% J.M. Schepers & M. van der Leest %%%
%% start script
clear all
clc
% close all

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

OBS_rank = rank(ctrb(LTI.A',LTI.C')); %Check observability
rank([eye(size(LTI.A,1))-LTI.A -LTI.B;LTI.C LTI.D]);

% disturbance matrixes
LTI.Cd=diag([0 0 0 0]);
LTI.Bd=[zeros(12,4)];
LTI.d=[0;0;0;0];  

%% initial values
Np = 10; % Prediction horizon
Nc = 1; % Control horiozn
simTime = 30;
T = simTime/T_pred; % sampling steps

x0= [pi/4 0 0 0 pi/8 0 0.2 0 1 0 0.3 0]';                    % inital conditions x(7) - z states

% Reference signal
PSI_ref = 0*ones(1,T);
Z_ref = 1*ones(1,T);
Z_ref = [0*ones(1,T/4) 1*ones(1,T/2) 0*ones(1,T/4)];
X_ref = 0*ones(1,T);
Y_ref = 0*ones(1,T);
r = [PSI_ref;Z_ref;X_ref;Y_ref];   % reference values psi(x5),z(x7),x(x9) and y (x11)

% QR weights
Q_weight = 10; % MPC Q state matrix
R_weight = 0.1;  % MPC R control input weight

Q_weight_LQ = 1; % weights LQ controller
R_weight_LQ =5; % weight input LQ controller

% d_dist = [zeros(1,T/10) 0.01*ones(1,3*T/10) zeros(1,6*T/10)]; %Dist on z
% %         zeros(1,T/10) 0.05*ones(1,3*T/10) zeros(1,6*T/10); %Dist on x
% %         zeros(1,T/10) 0.02*ones(1,3*T/10) zeros(1,6*T/10)]; %Dist on y

%% Initialize state and state estimation vectors
% stating matrices sizes
dim.nx = size(Ad,1); %number of states, 
dim.nu = size(Bd,2); %number of inputs, 
dim.ny = size(Dd,1);
dim.nd = 1; %number of disturbance inputs
dim.Np = Np; %prediction horizon
dim.Nc = Nc; %control horizon
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
t = zeros(1,T);                     % time vector
dtilde = 0.001;                    % disturbance initial value for augmented system
dtilde_ref = 0;


xaug = [x0;dtilde];

Vf = zeros(1,T);                % terminal cost sequence
l = zeros(1,T);                 % stage cost sequence


%% H matrix for MPC
% S=dare(LTI.A,LTI.B,Q,R);            % terminal cost

S = 1000*eye(dim.nx);              % beta method applying weight on terminal set
Qbar = blkdiag(kron(eye(dim.Np),Q),S);
Qbar_extnd = Qbar*kron(ones(dim.Np+1,1),eye(dim.nx));
Rbar = kron(R,eye(dim.Np));
Rbar_extnd = Rbar*kron(ones(dim.Np,1),eye(dim.nu));
Sbar = S;

[P,Z] = predmodgen1(LTI,dim);
             
H = (Z'*Qbar*Z + Rbar);

%% Augmented system for disturbance

Bdist = [0 0 0 0 0 0 1 0 0 0 0 0]';  % disturbance on Z state
Cdist = [0 1 0 0]';                 % output disturbance on Z state

A_aug = [LTI.A Bdist;zeros(dim.nd,dim.nx) eye(dim.nd)];

B_aug = [LTI.B;zeros(dim.nd,dim.nu)];
C_aug = [LTI.C Bdist];

HautusT = rank([eye(dim.nx)-Ad -Bdist;C_aug]); %Check controllability of augmented system


%% LQ Controller
Q_lq = Q_weight_LQ*eye(size(LTI.A));
R_lq = R_weight_LQ*eye(length(LTI.B(1,:)));
x_lq(:,1) = x0'; % Initial state LQ controller

[G_ref_lq,K_lq] = LQ_controller(G_ref,LTI,Q_lq,R_lq,T,T_pred);

%% MPC

u_limit = 4.9;
x_lim_vec = [0.5*pi 1000 0.5*pi 1000 pi 1000 1000 1000 1000 1000 1000 1000]';
x_lim_vec_full = repmat(x_lim_vec,[Np 1]);


for k = 1:T
    t(k) = (k-1)*T_pred;
    

    
    %%%%%%%% OTS %%%%%%%%%%
    Q_ots = eye(dim.nx);
    R_ots = eye(dim.nu);
    J_ots = blkdiag(Q_ots,R_ots);
    
    A_ots = [eye(dim.nx)-LTI.A -LTI.B; LTI.C zeros(dim.nx,dim.nu)]; 
    b_ots = [Bdist*dtilde_ref(:,k); G_ref*(y_ref(:,k) - Cdist*dtilde_ref(:,k))];
%     b_ots = [Bdist*dtilde(:,k); G_ref*(y_ref(:,k) - Cdist*dtilde(:,k))];
%     b_ots = [Bdist*d_dist(:,k); G_ref*(y_ref(:,k) - Cdist*d_dist(:,k))]; %     % this is old using changing disturbance

    
    opts = optimoptions('quadprog','Display','off');
    [xr_ur,~,exitflag] = quadprog(J_ots,zeros(dim.nx+dim.nu,1),[],[],A_ots,b_ots,[],[],[],opts);
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    x_ref_ots(:,k) = xr_ur(1:dim.nx);
    u_ref_ots(:,k) = xr_ur(1+dim.nx:dim.nx+dim.nu);

    x0 = xaug(1:dim.nx,k); %Estimated state
    f = x0'*P'*Qbar*Z - x_ref_ots(:,k)'*Qbar_extnd'*Z - u_ref_ots(:,k)'*Rbar_extnd';
 
    % compute control action
    cvx_begin quiet
        variable u_Np(4*Np)
        minimize ( (1/2)*quad_form(u_Np,H) + f*u_Np)
        %input constraints
        u_Np <=  u_limit*ones(4*Np,1);
        u_Np >= -u_limit*ones(4*Np,1);
        % state constraints
        Z(1:dim.nx*dim.Np,:)*u_Np <= -P(1:dim.nx*dim.Np,:)*x0 + x_lim_vec_full;
        Z(1:dim.nx*dim.Np,:)*u_Np >= -P(1:dim.nx*dim.Np,:)*x0 - x_lim_vec_full; 
    cvx_end
    
    u(:,k) = u_Np(1:dim.nu); % MPC control action
    
   % apply control action with disturbance
    xaug(:,k+1) = A_aug*xaug(:,k) + B_aug*u(:,k); % state / input disturbance
    yaug(:,k) = C_aug*xaug(:,k) + 0.01*randn(dim.ny,1);
    dtilde_ref(:,k+1) = dtilde_ref(:,k);
    dtilde(:,k+1) = xaug(13,k);
    
    %LQ controller
    u_lq(:,k) = -K_lq*x_lq(:,k);
    x_lq(:,k+1) = Ad*x_lq(:,k) + Bd*u_lq(:,k) + G_ref_lq*r(:,k) + 0.01*randn(dim.ny,1); % input disturbance
    y_lq(:,k) = Cd*x_lq(:,k); 
    
end

% states_trajectory: Nx12 matrix of trajectory of 12 states
states_trajectory_aug = yaug';
states_trajectory_lq = y_lq';
%%
figure
stairs(t,states_trajectory_aug(:,7))
hold on
stairs(t,states_trajectory_lq(:,7));
stairs(t,x_ref_ots(7,:),'k--')
xlabel('Time [s]')
ylabel('Altitude z[m]')
legend('Altitude with MPC - dist z[m]','Altitude with LQ - dist z[m]','reference altitude')

figure
plot(t,u(1,:))
hold on
plot(t,u_ref_ots(1,:))
title('Control input u1 vs reference')

% figure
% plot(t,dtilde_ref(1,1:end-1))
% hold on
% plot(t,d_dist)
% title('Estimated disturbance')
