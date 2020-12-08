clc;
clear;
load('matrix.mat');
%% 
y_sp=[0.4;0.8];
disturbance=[0.3;0.2];
rank([A B;C zeros(2,2)]);
A_bar=[A zeros(4,2);-C zeros(2,2)];
B_bar=[B;zeros(2,2)];
B_w_bar=[B;zeros(2,2)];
B_r_bar=[zeros(4,2);eye(2)];
C_bar=[C zeros(2,2)];
Q=[1 0 0 0 0 0;
   0 1 0 0 0 0;
   0 0 1 0 0 0;
   0 0 0 1 0 0;
   0 0 0 0 1 0;
   0 0 0 0 0 1]*10;
R=[1 0
   0 1]*0.5;
gamma=[A_bar -B_bar/R*B_bar'
    -Q -A_bar'];
[vector,value]=eig(gamma);
value=sum(value);
v=vector(:,find(real(value)<0));
P=v(7:12,:)/v(1:6,:);
K=real(inv(R)*B_bar'*P);
K1=K(:,1:4)
K2=K(:,5:6)
%% observe
syms l11 l12 l13 l14 l21 l22 l23 l24 s ;
dir_pole=[-12,-11.6,-10,-13];
polynomial=(s-dir_pole(1))*(s-dir_pole(2))*(s-dir_pole(3))*(s-dir_pole(4));
Ad_cof=double(coeffs(polynomial));
L=[l11 l12 l13 l14 ;l21 l22 l23 l24];
W=[C' A'*C' A'^2*C' A'^3*C'];
assert(rank(W(:,1:4))==4);
CC=W(:,1:4);
CC=[CC(:,1),CC(:,3),CC(:,2),CC(:,4)];
inv_C=inv(CC);
d1=2;
d2=2;
T=[inv_C(d1,:); 
    inv_C(d1,:)*A';
    inv_C(d1+d2,:);
    inv_C(d1+d2,:)*A'];
Aba=T*A'/(T);
Bba=T*C';
Aba(abs(Aba)<10^(-10))=0;
Bba(abs(Bba)<10^(-10))=0;
Am=Aba-Bba*L;
Ad=[ 0  1  0  0 ; 
     0  0  1  0 ;
     0  0  0  1 ;
    -Ad_cof(1:4) ];
 %% solve equetion
equation=Am==Ad;
L_num=solve(equation);
L_ans=struct2array(L_num);
L_ans=double(L_ans);
Lba=[L_ans(1:4);
    L_ans(5:8)];
L=Lba*T;
L=L'
L=real(L)