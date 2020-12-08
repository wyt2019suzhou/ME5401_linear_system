clc
clear
% syms s;
% lamda1=-1+0.5*i
% lamda2=-1-0.5*i
% polynomial=(s-lamda1)*(s-lamda2)
% Ad_cof=double(coeffs(polynomial));
% subplot(2,1,1)
% step(tf([1],fliplr([Ad_cof])))
% legend('pole=-1+-0.5i')
% lamda3=-1+10*i
% lamda4=-1-10*i
% polynomial2=(s-lamda3)*(s-lamda4)
% Ad_cof2=double(coeffs(polynomial2));
% subplot(2,1,2)
% step(tf([1],fliplr([Ad_cof2])));
% legend('pole=-1+-10i')
% clc;
% clear;
% load('matrix.mat')
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
lamda3=-1.28- 0.24i;
lamda4=-1.28+ 0.24i;
polynomial1=(s-lamda1)*(s-lamda2)*(s-lamda3)*(s-lamda4);
polynomial2=(s-lamda1)*(s-lamda2);
Ad_cof=double(coeffs(polynomial1));
Ad_cof2=double(coeffs(polynomial2));
step(tf([1],fliplr([Ad_cof])))
hold on 
step(tf([1],fliplr([Ad_cof2])))
legend('add extra pole','original system')

