function [G_ref_lq,K_lq,gain_LQ] = LQ_controller(G_ref,LTI,Q_lq,R_lq,T,T_pred)

G_ref_lq = G_ref;
x_lq = zeros(length(LTI.A(:,1)),T);    % state trajectory
u_lq = zeros(length(LTI.B(1,:)),T);    % control inputs
y_lq = zeros(length(LTI.C(:,1)),T);

% Compute feedback gain
[K_lq,~,~] = dlqr(LTI.A,LTI.B,Q_lq,R_lq,[]);
cl_lq = ss(LTI.A-LTI.B*K_lq,G_ref,LTI.C,[],T_pred);  % Computing state space
dcgain_cl_lq = dcgain(cl_lq);                        % Scale input

gain_LQ = dcgain_cl_lq;

G_ref_lq(5,1) = 1/dcgain_cl_lq(5,1);
G_ref_lq(7,2) = 1/dcgain_cl_lq(7,2);
G_ref_lq(9,3) = 1/dcgain_cl_lq(9,3);
G_ref_lq(11,4) = 1/dcgain_cl_lq(11,4);


end