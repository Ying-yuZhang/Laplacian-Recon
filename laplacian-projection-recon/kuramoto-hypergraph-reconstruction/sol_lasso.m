function A_i=sol_lasso(C,y,RHO,LAMBDA,t_max,MU)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%求解带有稀疏约束的线性方程组：min ||C*A_i-y||^2_2+LAMBDA|A_i|
%Copyright：马闯 chuang_m@126.com
%输入：C,y，其中CA_i=y
%      RHO,惩罚参数，一般取1
%      LAMBDA，一范数约束系数
%      t_max，最大迭代次数    
%      MU，   收敛误差      
%输出：A_i：方程组的解
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RHO=1;
% MU=0.001;
% LAMBDA=0.01;
% t_max=100000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[m,n]=size(C);%n表示未知量个数

%初始化
A_i=zeros(n,1);
V=zeros(n,1);
Z=zeros(n,1);
%事先处理
PHI=inv(C'*C+RHO*eye(n,n));
G=C'*y;
%迭代
for i=1:t_max
    A0=A_i;
    %更新A_i
    A_i=PHI*(G+RHO*(Z-V/RHO));
    %更新Z
    s1=(V/RHO+A_i)-LAMBDA/RHO;
    s2=-(V/RHO+A_i)-LAMBDA/RHO;
    Z=(s1+abs(s1))/2-(s2+abs(s2))/2;
    %更新V
    V=V+RHO*(A_i-Z);
    if max(abs(A_i-A0))<MU
       break;
    end
end
% A_i=Z;%把Z作为输出结果