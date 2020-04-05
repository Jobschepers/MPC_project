function [P,Z,W] =predmodgen1(LTI,dim)

%Prediction matrices generation
%This function computes the prediction matrices to be used in the
%optimization problem

%Prediction matrix from initial state
P=zeros(dim.nx*(dim.Np+1),dim.nx);
for k=0:dim.Np
    P(k*dim.nx+1:(k+1)*dim.nx,:)=LTI.A^k;
end

%Prediction matrix from input
Z=zeros(dim.nx*(dim.Np+1),dim.nu*(dim.Np));
for k=1:dim.Np
    for i=0:k-1
        Z(k*dim.nx+1:(k+1)*dim.nx,i*dim.nu+1:(i+1)*dim.nu)=LTI.A^(k-1-i)*LTI.B;
    end
end

W = zeros(dim.ny,dim.nu*dim.Np);
for i=0:dim.Np-1
    W(1:dim.ny,i*dim.nu+1:(i+1)*dim.nu) = LTI.C*LTI.A^(dim.Np-i-1)*LTI.B;
end
