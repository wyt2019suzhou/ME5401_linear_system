clc;
clear;
load('matrix.mat');
%% 
syms s 
for i=1:4
    if C(1,:)*A^(i-1)*B~=0
       degree1=i;
       break
    end
end
for i=1:4
    if C(2,:)*A^(i-1)*B~=0
       degree2=i;
       break
    end
end
Bstar=[C(1,:)*A^(degree1-1)*B
       C(2,:)*A^(degree2-1)*B];
F=inv(Bstar);
damp=0.8
wn=0.4
A_d=(s^2+2*damp*wn*s+wn^2);
A_d_den=double(fliplr(coeffs(A_d)));
% fA1=A+3*eye(4)
fA1=A+A_d_den(2)*eye(4)
K=F*[C(1,:)*fA1
     C(2,:)*fA1];
Bf=B*F;
Af = A-B*K
decouple_model=ss(Af,Bf,C,D)
% step(decouple_model)
p=pole(decouple_model)
H=C*inv(s*eye(4)-Af)*Bf
% t=0:0.1:30;
% x=zeros(4,1)
% u1=[ones(size(t,2),1),zeros(size(t,2),1)];
% lsim(decouple_model,u1,t,x);
% [vector,velue]=eig(Af)
% t=0:0.01:6;
% u0=zeros(size(t,2),2)
% [y,tout,x]=lsim(decouple_model,u0,t,x0);
% plot(t,x)