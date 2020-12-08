clc;
clear;
load('matrix.mat');
B1=[1 0;1 0;1 1;1 1];
xsp=[0;0.5;-0.4;0.3];
syms u1 u2 u3 u4 s
U=[u1;u2;u3;u4]
B2=[B B1];
eq1=inv(s*eye(4)-A)*B2*U;
eq2=subs(eq1,s,0);
ans=solve(eq2==xsp)
u1=double(ans.u1)
u2=double(ans.u2)
u3=double(ans.u3)
u4=double(ans.u4)
U1=[u1,u2]
U2=[u3,u4]