function [w,triangles]=generate_ER_HG(N,p,k2)
%底部是ER网络，三体是度不相关的均匀连接形成
%%%%%%输入%%%%%%%%%
% N 节点数; p 连边的概率；k2 每个节点的平均所属的三角形数；
%%%%%%输出%%%%%%%%%%%%
% w 邻接矩阵；triangles 3超边中的实三角形
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% N=10;
% p=0.5;
% k2=1;

%首先生成一个ER网
w=erdos_reyni(N,p);
w=spones(w);

%三体间连接概率与三体间度的相关性
p2=(2*k2)/((N-1)*(N-2));  
%生成三角形
% triangles_node=combnk(1:N,3);%三个点的所有可能排列组合
triangles_node = nchoosek(1:N,3);
triangles=triangles_node(rand(1,length(triangles_node))<p2,:);
end



