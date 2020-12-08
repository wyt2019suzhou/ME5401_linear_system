clc;
clear;
load('matrix.mat');
%% lqr
Q1=[1 0 0 0 
   0 1 0 0 
   0 0 1 0 
   0 0 0 1 ].*50;
%R=[1 0
   %0 1].*0.5;
Q2=[1 0 0 0 
   0 1 0 0 
   0 0 1 0 
   0 0 0 1 ];
R=[1 0
   0 1];
gamma1=[A -B/R*B'
    -Q1 -A'];
gamma2=[A -B/R*B'
    -Q2 -A'];
[vector1,value1]=eig(gamma1);
value1=sum(value1);
v1=vector1(:,find(real(value1)<0));
P1=v1(5:8,:)/v1(1:4,:);
K1=inv(R)*B'*P1;
[vector2,value2]=eig(gamma2);
value2=sum(value2);
v2=vector2(:,find(real(value2)<0));
P2=v2(5:8,:)/v2(1:4,:);
K2=inv(R)*B'*P2;
%% draw
t=0:0.01:5;
Af1=A-B*K1;
Af2=A-B*K2;
ss1=ss(Af1,B,C,D);
ss2=ss(Af2,B,C,D);
%x=zeros(4,1)
u0=zeros(size(t,2),2)
%u1=[ones(size(t,2),1),zeros(size(t,2),1)];
%u2=[zeros(size(t,2),1),ones(size(t,2),1)];
%[y1,tout1,x1]=lsim(ss1,u0,t,x0);
%plot(t,x1)
%legend('x1 Q=50I','x2 Q=50I','x3 Q=50I','x4 Q=50I')

%hold on
[y2,tout2,x2]=lsim(ss2,u0,t,x0);
plot(t,x2)
legend('x1 Q=I','x2 Q=I','x3 Q=I','x4 Q=I')
xlabel('time')
ylabel('state')
%figure(2)
%plot(t,y)
%legend('out1','out2')
%lsim(ss,u1,t,x);
%hold on
%lsim(ss1,u2,t,x);

%lsim(ss2,u2,t,x);
%legend('Q=10I','Q=I')