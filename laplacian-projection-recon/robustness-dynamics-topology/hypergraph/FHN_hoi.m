function dphi=FHN_hoi(t,phi,A,B,alpha,I)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%需要差分的函数
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = length(A); % 网络中振子的数量

u=phi(1:N,1);
v=phi((N+1):end,1);

% parameters of the units
a = 3;
b = 2;
c = 1;
epsilon = 1;
sigma = 1;

% ------------------------------
% 计算局部动力学项
% ------------------------------
f = c * (u - (u.^3)/3 - v);
g = c * (a*u - b*v);

% ------------------------------
% 计算网络耦合项
% ------------------------------
L=alpha(1)*A+alpha(2)*sum(B,3);
L=L-diag(sum(L,2));
coupling_u = epsilon * (L * u);
coupling_v = sigma * epsilon * (L * v);

% ------------------------------
% 合成导数
% ------------------------------
du = f + coupling_u+I;
dv = g + coupling_v;
dphi=[du;dv];

end

