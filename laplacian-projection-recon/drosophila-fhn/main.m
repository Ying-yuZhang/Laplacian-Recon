clear;clc;
load Net.mat
N = length(A);          %网络中振子的数量
M=N;
rng(42);

% initial conditions for the oscillators
xoold =10*rand(N,1);
yoold =10*rand(N,1);
X0=[xoold; yoold];

delta = 0.2;%驱动大小
[U, ~, V] = svd(rand(N)); % 进S行奇异值分解
di = diag(0.5 + rand(1, N)); % 生成非零奇异值，确保矩阵满秩
I = U * di * V'; % 重构满秩矩阵
% 获取矩阵 I 的最小值和最大值
I_min = min(I(:));
I_max = max(I(:));
% 将 I 的值映射到 [0, delta] 范围
I = 0 + (I - I_min) * (delta - 0) / (I_max - I_min);

tmax = 0.1;

T = [0 tmax];
% Accurate integration of the equation
options = odeset('abstol',1e-12,'reltol',1e-12);
[T,Xm]=ode45(@(t,x) FHN_hoi2(t,x,L,zeros(N,1)),T,X0,options);
X00=Xm(end,:);

f0=FHN_hoi2(0,X00(:),L,zeros(N,1));
X=X00;
theta=zeros(M,2*N);
D=zeros(M,2*N);
for j=1:M
    [T,Xm]=ode45(@(t,x) FHN_hoi2(t,x,L,I(:,j)),T,X0,options);
    X_m=Xm(end,:);
    f1=FHN_hoi2(0,X_m',L,I(:,j));
    X=[X;X_m];
    D(j,:)=(f1-f0-[I(:,j);zeros(N,1)])';
    theta(j,:)=X_m-X00;
end

if M<N
    X_score=sol_lasso(theta(:,1:N),D(:,1:N),1,1e-4,100000,10^-5);
else
    X_score=theta(:,1:N)\D(:,1:N);
end

ite_num=5;                  %迭代次数
cut=0.05;                   %邻接矩阵截断阈值
order=5;                    %最高阶数
delta=0.02;
init_alpha=Init_alpha_SC(X_score,order,cut); %估计的初始耦合强度

clique_sizes = 2:6;
true_data = {c2, c3, c4, c5, c6};  % 用评价性能
F1_Score=zeros(order,ite_num);

for ite=1:ite_num
    G = Extract_simplex(X_score,ite,init_alpha,delta,cut,order);

    [row, col] = find(tril(G.A_estimate, -1));  % 获取下三角部分的边
    rc2 = [col, row];  % 边列表

    pred_data = {rc2, G.A2_simplex_list, G.A3_simplex_list, G.A4_simplex_list, G.A5_simplex_list};

    for i = 1:length(clique_sizes)
        T = true_data{i};
        P = pred_data{i};
        F1_Score(i,ite) = calc_f1(P, T);
    end
end

save F1_Score.mat F1_Score
save G.mat G
