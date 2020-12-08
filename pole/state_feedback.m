clc;
clear;
load('matrix.mat');
%% reference model
syms s;
ts=20;%settling time
mp=0.1;%overshoot
pos_damp=log(1/mp)/sqrt(pi^2+(log(1/mp))^2);
pos=4/(ts);
damp=0.8;
wn=0.4;
lamda1=-damp*wn+wn*sqrt(1-damp^2)*i;
lamda2=-damp*wn-wn*sqrt(1-damp^2)*i;
lamda3=-1.5;
lamda4=-1.28;
polynomial=(s-lamda1)*(s-lamda2)*(s-lamda3)*(s-lamda4);
Ad_cof=double(coeffs(polynomial));
step(tf([1],fliplr([Ad_cof])));
%k = place(A,B,[lamda1 lamda2 lamda3 lamda4])
%% CCF
syms k11 k12 k13 k14 k21 k22 k23 k24 ;
K=[k11 k12 k13 k14; k21 k22 k23 k24];
W=[B A*B A^2*B A^3*B];
assert(rank(W(:,1:4))==4);
CC=W(:,1:4);
CC=[CC(:,1),CC(:,3),CC(:,2),CC(:,4)];
inv_C=inv(CC);
d1=2;
d2=2;
T=[inv_C(d1,:); 
    inv_C(d1,:)*A;
    inv_C(d1+d2,:);
    inv_C(d1+d2,:)*A];
Aba=T*A/(T);
Bba=T*B;
Aba(abs(Aba)<10^(-10))=0;
Bba(abs(Bba)<10^(-10))=0;
Am=Aba-Bba*K;
Ad=[ 0  1  0  0 ; 
     0  0  1  0 ;
     0  0  0  1 ;
     -Ad_cof(1:4)];
 %% solve equetion
equation=Am==Ad;
K_num=solve(equation);
K_ans=struct2array(K_num);
K_ans=double(K_ans);
Kba=[K_ans(1:4);
    K_ans(5:8)];
K1=Kba*T;
%% test
t=0:0.1:30;
Af=A-B*K1;
ss=ss(Af,B,C,D);
x=zeros(4,1)
%u0=zeros(size(t,2),2)
u1=[ones(size(t,2),1),zeros(size(t,2),1)];
%u2=[zeros(size(t,2),1),ones(size(t,2),1)];
%[y,tout,x]=lsim(ss,u0,t,x0);
%plot(t,x)
%legend('x1','x2','x3','x4')
%xlabel('time')
%ylabel('state')
%figure(2)
%plot(t,y)
%legend('out1','out2')
lsim(ss,u1,t,x);
%hold on
%lsim(ss,u2,t,x);
