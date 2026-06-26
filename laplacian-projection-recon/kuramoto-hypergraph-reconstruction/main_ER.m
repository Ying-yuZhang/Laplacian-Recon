function main_ER
N=20;   %节点数
k1=6;  %每个节点的平均所属的边数；
k2=2;   %每个节点的平均所属的三角形数；
% p=0.2449;  %底层图连边的概率；
p=k1/(N-1);
while true
[w,triangles]=generate_ER_HG(N,p,k2);
% full(w);
% triangles;
rk1=sum(sum(w))/N;
rk2=size(triangles,1)*3/N;
if round(rk1)==k1&&round(rk2)==k2
    break;
end
end


N=length(w);
A=full(w);
B=zeros(N,N,N);
for i=1:size(triangles,1)
    pe=perms(triangles(i,:));

    for idx = 1:size(pe, 1)
        B(pe(idx, 1), pe(idx, 2), pe(idx, 3))=1;
    end
end
eval(['save ERHG_',num2str(N),'_',num2str(round(k1)),'_',num2str(round(k2)),'_A.mat A'])
eval(['save ERHG_',num2str(N),'_',num2str(round(k1)),'_',num2str(round(k2)),'_B.mat B'])

