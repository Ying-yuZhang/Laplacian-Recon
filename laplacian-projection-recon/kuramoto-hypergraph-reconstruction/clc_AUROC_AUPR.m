function [AUROC, AUPR]=clc_AUROC_AUPR(x,x0)
%     x0=[1,0,1,0,0,1]
%     x=[0.5,0.1,0,0,0,0]      
    P=sum(x0==1);            %1的个数
    N=sum(x0==0);            %0的个数
    y=sort(unique(x),'descend');  %去掉重复的，并从大到小排序
    n1=length(y);
    TPR=zeros(1,n1);
    FPR=zeros(1,n1);
    Precision=zeros(1,n1);
    for i=1:n1
      id=find(x>=y(i));
      TPR(i)=sum(x0(id)==1)/P;
      FPR(i)=sum(x0(id)==0)/N;
      Precision(i)=sum(x0(id)==1)/length(id);
    end
    
   Recall=TPR;
   FPR=[0,FPR,1]; 
   TPR=[0,TPR,1];
   Recall=[0,Recall,1];
   Precision=[Precision(1),Precision,Precision(end)];
   AUROC=trapz(FPR,TPR);
   AUPR=trapz(Recall,Precision);
