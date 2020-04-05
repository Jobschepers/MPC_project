%% STABILITY OF LINEARIZED QUADROTOR AND INVERSE PENDULUM SYSTEM
 
% Given a set of linearized system matrices, this code determines the
% terminal set X_f which guarantees asymptotic stability. Finally, the
% value of X_N (calygraph X) can be be determined for a given prediction
% horizon N

% This script is dowloaded from Alex Berndt on GitHub https://github.com/alexberndt/Quadrotor-Balancing-Pendulum-MPC
% and is modified for the project of MPC SC42125 - Tu Delft
% by Job Schepers and Max van der Leest for the given LTI system
% Linearizedmodel.mat 

% 

%%
clc
clear
addpath('functions/');

disp('------------------------------------------------------------------');
disp('          STABILITY ANALYSIS OF FLYING QTW UAV ');
disp('');
disp('------------------------------------------------------------------');

%% INITIALIZATION
fprintf('\tinitializing ... \n');

% DISCRETIZE SYSTEM
load('linearizedmodel.mat')


T_pred = 0.1;% sample time
lambda = 0; % time delay
[Ad,Bd,Cd,Dd] = c2dt(sys.A,sys.B,sys.C,T_pred,lambda);
LTI.A= Ad;
LTI.B= Bd;
LTI.C= Cd;
LTI.D= Dd;
fprintf('\tdone!\n');

%% DETERMINE OPTIMAL LQR GAIN FOR MPC COST FUNCTION
fprintf('\tdetermining LQR optimal control action ... \n');

Q = 10*eye(size(Ad));           % state quadratic cost 
R = 0.1*eye(length(Bd(1,:)));   % input quadratic cost

[X,L,K_LQR] = dare(Ad,Bd,Q,R);  % determine LQR gain for unconstrained system

A_K = Ad-Bd*K_LQR;              % closed-loop LQR system
eigvals_A_K = eig(A_K);         % determine closed-loop eigenvalues

fprintf('\tdone!\n');
%% DETERMINE INVARIANT ADMISSIBLE SET X_f
fprintf('\testimating X_f invariant set \n');

dim.nx = 12;        %
dim.nu = 4;         %

u_limit = 4.9;      % bound on control inputs
x_limit = 0.5*pi;        % bound on states
x_lim_vec = [0.5*pi 1000 0.5*pi 1000 pi 1000 1000 1000 1000 1000 1000 1000]';        % bound on states

fprintf('\t - defining state and input constraints \n');
Fu = kron(eye(dim.nu),[1; -1]);
Fx = kron(eye(dim.nx),[1; -1]);

fu = u_limit*ones(2*dim.nu,1);
fx = x_limit*ones(2*dim.nx,1);

f = [fu; fx];
F = blkdiag(Fu,Fx);

s = size(F,1);

C_aug = [K_LQR; eye(dim.nx)];


%% Construct terminal set
[Xf_set_H,Xf_set_h] = max_output_set(A_K,K_LQR,u_limit,x_lim_vec);

fprintf('\tSuccesfully constructed terminal set X_f\n');

%% CHECK CONSTRAINT

% given X_f calculated before, we can test if a given state x0 is within
% this set or not. This is done here:

% state x0
%       <-- state names
x0= [0 0 0 0 0 0 0.8 0 0.1 0 0 0]'; % <-- state values

% check if the x-location is within X_f
inSet = all(Xf_set_H*x0 <= Xf_set_h)

%% TEST IF MPC LAW GUIDES TOWARDS X_f IN N-STEPS

% Having constructed X_f, we empirically construct a state X_N (calygraph X) 
% to determine which states can be steered to the terminal set X_f in N 
% steps.

fprintf('\tDetermining X_N emperically for different Beta values ...\n');

N = 20;                     % prediction horizon
fprintf('\t - N = %i (prediction horizon) \n',N);

x = zeros(dim.nx,N+1);      % state trajectory

% initial state
x(:,1) = x0;

% tuning weights
Q = 10*eye(dim.nx);         % state cost
R = 0.1*eye(dim.nu);        % input cost

% terminal cost = unconstrained optimal cost (Lec 5 pg 6)
[S,~,~] = dare(Ad,Bd,Q,R);  % terminal cost 
S = 1000*eye(dim.nx);

% determine prediction matrices
% Qbar = kron(Q,eye(N));
Qbar = blkdiag(kron(eye(N),Q),S);
Rbar = kron(R,eye(N));
Sbar = S;

LTI.A = Ad;
LTI.B = Bd;
LTI.C = Cd;

dim.N = N;
dim.Np = N;
dim.nx = size(Ad,1);
dim.nu = size(Bd,2);
dim.ny = size(Cd,1);

% [P,Z,W] = predmodgen_exp(LTI,dim);
[P,Z,W] = predmodgen1(LTI,dim); % incorrect dimensions matrixes

             
% define input constraints
lb = -u_limit*ones(4*N,1);
ub = u_limit*ones(4*N,1);

Aeq = Xf_set_H;
beq = Xf_set_h;

% Define IC_test_vals as the set of initial conditions to consider
r = 1;
s = 1;

res = 10;
bnd_z = 10; % ~~~~~~~~~ Changed from 0.2
bnd_psi = 2*pi; % ~~~~~~~~~ Changed from 0.2

mat = zeros(res,res,3);
% matPedro = zeros(res,res);
beta_i = 1;

for betaVal = [0.1 1 10]
    
    Sbar = betaVal*S;
    fprintf('\t - Beta = %d \n',betaVal)

    H = (Z'*Qbar*Z + Rbar + 2*W'*Sbar*W);
    d = (x0'*P'*Qbar*Z + 2*x0'*(Ad^N)'*Sbar*W)';
    
    r = 1;

    % loop through different initial conditions 
    for initial_psi = linspace(-bnd_psi,bnd_psi,res)
        s = 1;
        for initial_z = linspace(-bnd_z,bnd_z,res)
            x(:,1) = [0 0 0 0 0 initial_psi 0 initial_z 0 0 0 0]';
%             x(:,1) = [0 0 0 0 0 0 0 initial_z 0 0 initial_psi 0]';
            for k = 1:N     
                % define current state position
                x_current = x(:,k);
                d = (x_current'*P'*Qbar*Z + 2*x_current'*(Ad^N)'*Sbar*W)';
                warning off
                % solve QP problem
                opts = optimoptions('quadprog','Display','off');
                uopt = quadprog(H,d,[],[],[],[],lb,ub,[],opts);
                warning on

                % obtain optimal control action at k=0
                u(:,k) = uopt(1:4);

                % apply control action
                x(:,k+1) = Ad*x(:,k) + Bd*u(:,k);
            end

            %disp('Does x(N) belong to X_f?:');
            %disp(initial_r);
            %disp(initial_s);
            inSet = all(Xf_set_H*x(:,N+1) <= Xf_set_h);
            %disp(inSet);

            mat(r,s,beta_i) = inSet;
            % matPedro(r,s) = all(H_pedro*x(:,N+1) <= h_pedro);

            s = s+1;
        end
        r = r+1;
    end
    
    beta_i = beta_i + 1;
end

fprintf('\tdone!\n');

% plot the initial conditions which were steered towards the terminal set
% X_f within N steps
%%
x_as = linspace(-bnd_z,bnd_z,res);
y_as = linspace(-bnd_psi,bnd_psi,res);
c = mat(:,:,2);
mymap = [1 0 0;
         0 1 0];
figure;
title('Emperical plot of X_N');
imagesc(x_as,y_as,mat(:,:,2));
colormap(mymap)
%  axis([x_as(1) x_as(end) y_as(1) y_as(end)])
xlabel('Initial values Z[m]')
ylabel('Initial values $\psi$ [rad]','interpreter','latex');

%%

fprintf('\tplotting results ... \n');

figure;
title('Emperical plot of X_N');
subplot(3,1,1);
imagesc(x_as,y_as,mat(:,:,1));
colormap(mymap)
xlabel('Initial values Z[m]')
ylabel({'$\beta$ = 0.1*s',' Initial values $\psi$ [rad]'},'interpreter','latex');

subplot(3,1,2);
imagesc(x_as,y_as,mat(:,:,2));
colormap(mymap)
xlabel('Initial values Z[m]')
ylabel({'$\beta$ = 1*s',' Initial values $\psi$ [rad]'},'interpreter','latex');

subplot(3,1,3);
imagesc(x_as,y_as,mat(:,:,3));
colormap(mymap)
xlabel('Initial values Z[m]')
ylabel({'$\beta$ = 10*s ',' Initial values $\psi$ [rad]'},'interpreter','latex');

fprintf('\tdone! \n');

disp('------------------------------------------------------------------');
disp('            INVARIANT CONTROL ADMISSIBLE SET X_f ');
disp('');
disp('------------------------------------------------------------------');
disp('Hx < h defines the set X_f');
disp('See workspace for H and h as Xf_set_H and Xf_set_h respectively');

%% PLOT THE CONSTRAINT SET X_f

% H_xyz = Xf_set_H(:,[3,9,13]);
% 
% 
% figure(3);
% plot3(H_xyz(:,1),H_xyz(:,2),H_xyz(:,3),'*m');
% grid();
% 
% 
% [K,V] = convhull(H_xyz(:,1),H_xyz(:,2))
% hold on
% plot(H_xyz(K,1),H_xyz(K,2),'r-');


% xx = -1:.05:1;
% yy = abs(sqrt(xx));
% [x,y] = pol2cart(xx,yy);
% k = convhull(x,y,x+1);
% figure(9);
% plot(x(k),y(k),'r-',x,y,'b*')

%% PLOT OF REPORT
% %% TEST IF MPC LAW GUIDES TOWARDS X_f IN N-STEPS
% 
% % Having constructed X_f, we empirically construct a state X_N (calygraph X) 
% % to determine which states can be steered to the terminal set X_f in N 
% % steps.
% 
% fprintf('\tDetermining X_N emperically ...\n');
% 
% N = 10;                     % prediction horizon
% fprintf('\t - N = %i (prediction horizon) \n',N);
% 
% x = zeros(dim.nx,N+1);      % state trajectory
% 
% % initial state
% x0 = [0.05 0 0.1 0 0 0  0.05 0 0.4 0 0 0  0.2 0  0.3 0]';
% x(:,1) = x0;
% 
% % tuning weights
% Q = 10*eye(dim.nx);         % state cost
% R = 0.1*eye(dim.nu);        % input cost
% 
% % terminal cost = unconstrained optimal cost (Lec 5 pg 6)
% [S,~,~] = dare(Ad,Bd,Q,R);  % terminal cost % OLD: S = 10*eye(size(A));
% 
% % determine prediction matrices
% Qbar = kron(Q,eye(N));
% Rbar = kron(R,eye(N));
% Sbar = S;
% 
% LTI.A = Ad;
% LTI.B = Bd;
% LTI.C = Cd;
% 
% dim.N = N;
% dim.nx = size(Ad,1);
% dim.nu = size(Bd,2);
% dim.ny = size(Cd,1);
% 
% [P,Z,W] = predmodgen(LTI,dim);
%              
% % define input constraints
% u_limit = 0.1;
% lb = -u_limit*ones(4*N,1);
% ub = u_limit*ones(4*N,1);
% 
% Aeq = Xf_set_H;
% beq = Xf_set_h;
% 
% % Define IC_test_vals as the set of initial conditions to consider
% r = 1;
% s = 1;
% 
% res = 10;
% bnd_x = 0.8;
% bnd_y = 0.8;
% 
% mat_plot = zeros(res,res);
% betaVal = 1;
%     
% Sbar = betaVal*S;
% 
% H = (Z'*Qbar*Z + Rbar + 2*W'*Sbar*W);
% d = (x0'*P'*Qbar*Z + 2*x0'*(Ad^N)'*Sbar*W)';
% 
% r = 1;
% 
% % loop through different initial conditions 
% for initial_r = linspace(-bnd_x,bnd_x,res)
%     s = 1;
%     for initial_s = linspace(-bnd_y,bnd_y,res)
%         x(:,1) = [initial_r initial_s 0 0 0 0  0 0 0 0 0 0  0 0  0 0]';
%         for k = 1:N     
%             % define current state position
%             x_current = x(:,k);
%             d = (x_current'*P'*Qbar*Z + 2*x_current'*(Ad^N)'*Sbar*W)';
%             warning off
%             % solve QP problem
%             opts = optimoptions('quadprog','Display','off');
%             uopt = quadprog(H,d,[],[],[],[],lb,ub,[],opts);
%             warning on
% 
%             % obtain optimal control action at k=0
%             u(:,k) = uopt(1:4);
% 
%             % apply control action
%             x(:,k+1) = Ad*x(:,k) + Bd*u(:,k);
%         end
%         inSet = all(Xf_set_H*x(:,N+1) <= Xf_set_h);
%         mat_plot(r,s) = inSet;
%         % matPedro(r,s) = all(H_pedro*x(:,N+1) <= h_pedro);
%         s = s+1;
%     end
%     r = r+1;
% end
% 
% fprintf('\tdone!\n');
% 
% % plot the initial conditions which were steered towards the terminal set
% % X_f within N steps
% 
% fprintf('\tplotting results ... \n');
% figure(2);
% hold on;
% 
% r = 1;
% for r_plot = linspace(-bnd_x,bnd_x,res)
%     s = 1;
%     for s_plot = linspace(-bnd_y,bnd_y,res)
%         if (mat_plot(r,s) == 1)
%             plot(r_plot,s_plot,'sg','MarkerFaceColor','g','MarkerEdgeColor','g','MarkerSize',8);
%         else
%             plot(r_plot,s_plot,'sr','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',8);
%         end
%         s = s+1;
%     end
%     r = r+1;
% end
% xlabel('r');
% ylabel('$\dot{r}$','interpreter','latex');
% grid on;
% xlim([-1.05*bnd_x 1.05*bnd_x]);
% ylim([-1.05*bnd_y 1.05*bnd_y]);
% fprintf('\tdone! \n');







