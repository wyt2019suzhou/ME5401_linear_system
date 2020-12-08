clc;
clear;
load('matrix.mat');
system=ss(A,B,C,D)
G=tf(system)
syms s
N=[-0.7587*s^3 + 94.85*s^2 + 1251*s + 2815,-0.01574*s^3 + 36.07*s^2 + 465.8*s + 1060;-0.02941*s^3 + 34.36*s^2 + 382.8*s + 817.6,0.03102*s^3 + 13.1*s^2 + 162.3*s + 556.4]    
Kd=inv(N)*det(N)
ds=s^4 + 36.42*s^3 + 548.8*s^2 + 3381*s + 6013
res=(det(N)/ds)*eye(2) %GKd

thisnow=res;
answer=gcd(thisnow(1,1),thisnow(2,2));
Ks=[1/answer,0;0,1/answer]
K=Kd*Ks;
H=inv(eye(2)+ res*Ks)*res*Ks
s = 0
subs(s)
K1=eval(K)
H1=eval(H)
% syms s
% N=[(s+11.06)*(s+1.625)*(s^2 + 23.22*s + 324.9),-0.0068629*(s-5806)*(s+12.52)*(s+1.758);
%     -0.046374 *(s-716.5) *(s+9.127) *(s+3.743),(s+10.07) *(s+1.906) *(s^2 + 24.36*s + 241)]
% % G=N/((s+10.25)* (s+1.599) *(s^2 + 24.47*s + 231))
% G=[((s+11.06) *(s+1.625) *(s^2 + 23.22*s + 324.9))/((s+10.25)* (s+1.599)* (s^2 + 24.47*s + 231)),(-0.0068629* (s-5806)* (s+12.52)* (s+1.758))/((s+10.25)* (s+1.599)* (s^2 + 24.47*s + 231));
%     (-0.046374 *(s-716.5) *(s+9.127) *(s+3.743))/((s+10.25)* (s+1.599)* (s^2 + 24.47*s + 231)),(s+10.07) *(s+1.906) *(s^2 + 24.36*s + 241)/ ((s+10.25)* (s+1.599)* (s^2 + 24.47*s + 231))]
% Kd=inv(N)*det(N)
% % Kd=simplify(Kd)
% ds=(s+10.25) *(s+1.599) *(s^2 + 24.47*s + 231)
% res=(det(N)/ds)*(eye(2))^2 %GKd
% thisnow=res;
% answer=gcd(thisnow(1,1),thisnow(2,2));
% Ks=[1/answer,0;0,1/answer]
% K=Kd*Ks;
% % H=inv(eye(2)+ G*K)*G*K
% H=inv(eye(2)+ res*Ks)*res*Ks
% s = 0
% subs(s)
% K1=eval(K)
% H1=eval(H)