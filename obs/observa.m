clc;
clear;
load('matrix.mat');
%% LQR
syms s 
Q=[1 0 0 0 
   0 1 0 0 
   0 0 1 0 
   0 0 0 1 ].*10;
R=[1 0
   0 1].*0.5;
gamma=[A -B/R*B'
    -Q -A'];
[vector,value]=eig(gamma);
value=sum(value);
v=vector(:,find(real(value)<0));
P=v(5:8,:)/v(1:4,:);
K1=inv(R)*B'*P;
K1=real(K1)
Af=A-B*K1;
ss=ss(Af,B,C,D);
orig_pole=pole(ss);
dir_pole=[-9,-9,-9,-9];
polynomial=(s-dir_pole(1))*(s-dir_pole(2))*(s-dir_pole(3))*(s-dir_pole(4));
Ad_cof=double(coeffs(polynomial));
%% full order Observer LQR method
Q1=[1 0 0 0 
   0 1 0 0 
   0 0 1 0 
   0 0 0 1 ].*3.4;
R1=[1 0 
    0 1 ]*0.5;

phi1 = [A' -C'/R1*C;
            -Q1 -A];
[vector,value]=eig(phi1);
value=sum(value);
v=real(vector(:,find(real(value)<0)));
P=v(5:8,:)/v(1:4,:);
L1=inv(R1)*C*P;
L1=L1';
L1=real(L1)
A1=A-L1*C;
[vector,value]=eig(A1)
%% full order pole placement
syms l11 l12 l13 l14 l21 l22 l23 l24  ;
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
L2=Lba*T;
L2=L2'
L2=real(L2)
% %% reduce order
% syms t13 t14 t23 t24 g11 g12 g21 g22 
% Ob=[C;C*A;C*A*A;C*A*A*A];
% Obr=rank(Ob);
% D1=-1.*eye(2);
% G=[g11 g12 
%    g21 g22];
% T=[1 0 t13 t14;1 1 t23 t24]
% eq=D1*T-T*A+G*C==0;
% answer=solve(eq);
% answer=double(struct2array(answer));
% G=[answer(1:2)
%     answer(3:4)];
% T=[1 0 answer(5) answer(6);1 1 answer(7) answer(8)];
% E=T*B;
% H=inv([C;T]);