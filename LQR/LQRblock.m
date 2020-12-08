clc;
clear;
load('matrix.mat');
%% lqr
Q=[1 0 0 0 
   0 1 0 0 
   0 0 1 0 
   0 0 0 1 ]*100;
R=[1 0
   0 1];
%[K1,~,P]=lqr(A,B,Q,R)
gamma=[A -B/R*B'
    -Q -A'];
[vector,value]=eig(gamma);
value=sum(value);
v=vector(:,find(real(value)<0));
P=v(5:8,:)/v(1:4,:);
K1=real(inv(R)*B'*P);
%% draw
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