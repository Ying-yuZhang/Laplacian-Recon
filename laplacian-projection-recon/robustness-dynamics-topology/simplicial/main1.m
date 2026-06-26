name=char('ERSC_100_8_2','SWSC_100_8_2','SFSC_100_8_2','colony2','hypertext2009');
idlen=size(name,1);
perfor=zeros(9,5,idlen,5);
rng(42);
for l=1:5
    for id=1:idlen
        eval(['load ',strtrim(name(id,:)),'_A_',num2str(l),'.mat A'])% 邻接矩阵
        eval(['load ',strtrim(name(id,:)),'_B_',num2str(l),'.mat B'])% 高阶交互，2-单纯形

        N = length(A);          %网络中振子的数量
        %number of unknowns for each node
        H=(N-1)+(N-1)*(N-2)/2;
        omega = ones(N,1);      %固有频率
        alpha = [0.3,0.7];
        

        % initial conditions for the oscillators
        X0 = 0.5*rand(N,1);

        delta = 0.2;%驱动大小
        [U, ~, V] = svd(rand(N)); % 进S行奇异值分解
        di = diag(0.5 + rand(1, N)); % 生成非零奇异值，确保矩阵满秩
        I = U * di * V'; % 重构满秩矩阵
        % 获取矩阵 I 的最小值和最大值
        I_min = min(I(:));
        I_max = max(I(:));
        % 将 I 的值映射到 [0, delta] 范围
        I = 0 + (I - I_min) * (delta - 0) / (I_max - I_min);
        M = N;

        dx=zeros(M,N);

        tmax = 1;

        T = [0 tmax];
        % Accurate integration of the equation
        options = odeset('abstol',1e-12,'reltol',1e-12);
        [EdgeList, TriangleList] = convert_to_list(A, B);
        EdgeList0=[1 1; 2 2];
        TriangleList0=[1 1 1; 2 2 2];
        [T,Xm]=ode45(@(t,x) kuramoto_hoi2(t,x,EdgeList,TriangleList,omega,alpha,0),T,X0,options);
        X00=Xm(end,:);
        f0=kuramoto_hoi2(0,X00(:),EdgeList,TriangleList,omega,alpha,0);
        f=f0';
        X=X00;
        theta=zeros(M,N);
        for j=1:M
            [T,Xm]=ode45(@(t,x) kuramoto_hoi2(t,x,EdgeList,TriangleList,omega,alpha,I(:,j)),T,X0,options);
            X_m=Xm(end,:);
            f1=kuramoto_hoi2(0,X_m',EdgeList,TriangleList,omega,alpha,I(:,j));
            f=[f;(f1-I(:,j))'];
            X=[X;X_m];
            dx(j,:)=(f1-f0-I(:,j))';
            theta(j,:)=X_m(1:N)-X00;
        end
        D=dx(:,1:N);

        if M<N
            X_score=sol_lasso(theta,D,1,1e-4,100000,10^-5);
        else
            X_score=theta\D;
        end


        ite_num=10;                  %迭代次数
        cut=0.1;                   %邻接矩阵截断阈值
        order=2;                    %最高阶数
        delta=0.02;
        init_alpha=Init_alpha_SC(X_score,order,cut); %估计的初始耦合强度

        G = Extract_simplex(X_score,ite_num,init_alpha,delta,cut,order);

        %对比算法this
        XX=X(:,1:N)';
        YY=f(:,1:N)';
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


        % Zero order method
        Phi = zeros(nt,H);
        ff=f-omega';

        for i=1:N
            %pairwise interactions j=1:N-1
            Phi(:,1:N-1) = k*(sin(X(I1,[1:i-1 i+1:N])-X(I1,i)));
            AA(1:N-1,i)=A(i,[1:i-1 i+1:N]);
            %h.o.i. terms j=N:H
            vtemp=zeros(nt,(N-1)*(N-2)/2);
            itemp=1;
            for ii1=[1:i-1 i+1:N]
                for jj1=[1:i-1 i+1:N]
                    if jj1>ii1
                        vtemp(:,itemp)=kD*(sin(X(I1,ii1)+X(I1,jj1)-2*X(I1,i)));
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


            Yi = ff(I1,i);
            Zi0 = lsqminnorm(Phi,Yi,1e-12);
            Zi2 = lsqnonneg(Phi,Yi);

            %%%%%signal_lasso
            w0 = zeros(H,1);
            alpha1=0.001;
            alpha2=0.001;
            max_iters_ = 50; %50 for a quick check, 50000 for accurate results
            intercept_=0;
            Zi3 = signal_lasso(Phi,Yi,alpha1,alpha2,max_iters_,w0,intercept_);

            A00(:,i) = Zi0;
            A2(:,i) = Zi2;
            A3(:,i) = Zi3;
        end
        ours_AA=convert_to_AA(G.A_estimate, G.A2_estimate);
        eval(['save ',strtrim(name(id,:)),'_IKuramoto_',num2str(l),'.mat ours_AA A00 A2 A3 this_AA AA'])
        perfor(:, :, id, l) = clc_alg_performance(AA, A00, A2, A3, this_AA, ours_AA); 
        
    end
end
save IKuramoto_SC_perfor.mat perfor
