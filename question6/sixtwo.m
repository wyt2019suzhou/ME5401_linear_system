clc;
clear;
load('matrix.mat');
syms s u1 u2
U=[u1;u2]
a=5
b=4
c=8
d=0
w=[a+1 0 0 0;
    0 b+1 0 0;
    0 0 c+1 0;
    0 0 0 d+1]
xsp=[0;0.5;-0.4;0.3];
xs=inv(s*eye(4)-A)*B*U;
xs2=subs(xs,s,0);
J=1/2*(xs2-xsp)'*w*(xs2-xsp)
jacob = jacobian(J, [u1 u2])
ans=solve(jacob==0)
u1=double(ans.u1)
u2=double(ans.u2)
u=[u1;u2]
