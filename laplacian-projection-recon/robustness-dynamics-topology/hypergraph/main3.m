name=char('ERHG_100_6_2','SWHG_100_6_2','SFHG_100_6_2','high_school','InVS13');
idlen=size(name,1);
perfor=zeros(9,5,idlen,5);
rng(42);

for l=1:5
    for id=1:idlen
        eval(['load ',strtrim(name(id,:)),'_A_',num2str(l),'.mat A'])% 邻接矩阵
        eval(['load ',strtrim(name(id,:)),'_B_',num2str(l),'.mat B'])% 高阶交互，2-单纯形
        N = length(A);              %网络中振子的数量
        %number of unknowns for each node
        H=(N-1)+(N-1)*(N-2)/2;
        alpha = [0.3,0.7];
        M=N;


        % initial conditions for the oscillators
        xoold =rand(N,1)+10;
        yoold =rand(N,1)+10;
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

        tmax = 1;

        T = [0 tmax];
        % Accurate integration of the equation
        options = odeset('abstol',1e-12,'reltol',1e-12);
        [T,Xm]=ode45(@(t,x) FHN_hoi(t,x,A,B,alpha,zeros(N,1)),T,X0,options);
        X00=Xm(end,:);
        f0=FHN_hoi(0,X00(:),A,B,alpha,zeros(N,1));
        X=X00;
        theta=zeros(M,N);
        D=zeros(M,N);
        for j=1:M
            [T,Xm]=ode45(@(t,x) FHN_hoi(t,x,A,B,alpha,I(:,j)),T,X0,options);
            X_m=Xm(end,:);
            f1=FHN_hoi(0,X_m',A,B,alpha,I(:,j));
            X=[X;X_m];
            D(j,:)=(f1(1:N)-f0(1:N)-I(:,j))';
            theta(j,:)=X_m(1:N)-X00(1:N);
        end

        if M<N
            X_score=sol_lasso(theta,D,1,1e-4,100000,10^-5);
        else
            X_score=theta\D;
        end

        ite_num=10;                  %迭代次数
        cut=0.1;                   %邻接矩阵截断阈值
        order=2;                    %最高阶数
        delta=0.02;
        init_alpha=Init_alpha(X_score,order,cut);%[0.3 0.7];
        G = Extract_hypergraph(X_score,ite_num,init_alpha,delta,cut,order);

        %对比算法this
        this_f=zeros(size(X));
        for j=1:size(X,1)
            f1 = FHN_hoi(0,X(j,:)',A,B,alpha,zeros(N,1));
            this_f(j,:) = f1';
        end
        XX=X(:,1:N)';
        YY=this_f(:,1:N)';
        save('this_input.mat', 'XX', 'YY');
        status = system('julia this_interface.jl');

        if status ~= 0
            error('Julia脚本运行失败，请检查路径或代码');
        end

        load this_result.mat
        this_AA = This_convert_to_AA(this_EdgeList, this_TriangleList,N);

        %对比算法NNLS(耦合函数已知)
        k=0.3;
        kD=0.7;

        I1 = 2:size(X,1);     % Internal steps
        nt = size(I1,2);
        A00=zeros(H,N); %reconstructed matrix with OLS
        A2=zeros(H,N); %reconstructed matrix with NNLS
        A3=zeros(H,N); %reconstructed matrix with signal lasso
        AA=zeros(H,N); %"true" matrix

        [EdgeList, TriangleList] = convert_to_list(A, B);
        EdgeList0=zeros(size(A));
        TriangleList0=zeros(size(B));
        f=zeros(size(X));

        for j=1:size(X,1)
            f1 = FHN_hoi(0,X(j,:)',A,B,alpha,zeros(N,1))-FHN_hoi(0,X(j,:)',EdgeList0,TriangleList0,alpha,zeros(N,1));
            f(j,:) = f1';
        end

        % Zero order method
        Phi = zeros(nt,H);

        for i=1:N
            %pairwise interactions j=1:N-1
            Phi(:,1:N-1) = k*(X(I1,[1:i-1 i+1:N])-X(I1,i));
            AA(1:N-1,i)=A(i,[1:i-1 i+1:N]);
            %h.o.i. terms j=N:H
            vtemp=zeros(nt,(N-1)*(N-2)/2);
            itemp=1;
            for ii1=[1:i-1 i+1:N]
                for jj1=[1:i-1 i+1:N]
                    if jj1>ii1
                        vtemp(:,itemp)=kD*(X(I1,ii1)+X(I1,jj1)-2*X(I1,i));
                        RowIdx = find(ismember(TriangleList, [i ii1 jj1],'rows'));
                        if (~isempty(RowIdx))
                            AA(N-1+itemp,i)=1;
                        end
                        RowIdx = find(ismember(TriangleList, [ii1 i jj1],'rows'));
                        if (~isempty(RowIdx))
                            AA(N-1+itemp,i)=1;
                        end
                        RowIdx = find(ismember(TriangleList, [ii1 jj1 i],'rows'));
                        if (~isempty(RowIdx))
                            AA(N-1+itemp,i)=1;
                        end
                        itemp=itemp+1;
                    end
                end
            end
            Phi(:,N:N-1+(N-1)*(N-2)/2)=vtemp;

            Yi = f(I1,i);
            Zi0 = lsqminnorm(Phi,Yi,1e-12);
            Zi2 = lsqnonneg(Phi,Yi);

            %%%%%signal_lasso
            w0 = zeros(H,1);
            alpha1=0.001;
            alpha2=0.001;
            max_iters_ = 50;%50000; %50 for a quick check, 50000 for accurate results
            intercept_=0;
            Zi3 = signal_lasso(Phi,Yi,alpha1,alpha2,max_iters_,w0,intercept_);

            A00(:,i) = Zi0;
            A2(:,i) = Zi2;
            A3(:,i) = Zi3;
        end
        ours_AA=convert_to_AA(G.A_estimate, G.A2_estimate);
        eval(['save ',strtrim(name(id,:)),'_FHN_',num2str(l),'.mat ours_AA A00 A2 A3 this_AA AA'])

        perfor(:, :, id, l) = clc_alg_performance(AA, A00, A2, A3, this_AA, ours_AA); 
    end
end
save FHN_perfor.mat perfor
