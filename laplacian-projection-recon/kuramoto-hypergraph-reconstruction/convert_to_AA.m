function AA = convert_to_AA(A, B)
N = length(A);              %网络中振子的数量
%number of unknowns for each node
H=(N-1)+(N-1)*(N-2)/2;
AA=zeros(H,N); %"true" matrix

for i=1:N
    %pairwise interactions j=1:N-1
    AA(1:N-1,i)=A(i,[1:i-1 i+1:N]);
    %h.o.i. terms j=N:H
    itemp=1;
    for ii1=[1:i-1 i+1:N]
        for jj1=[1:i-1 i+1:N]
            if jj1>ii1
                if (B(i ,ii1, jj1)==1)
                    AA(N-1+itemp,i)=1;
                end
                if (B(ii1 ,i ,jj1)==1)
                    AA(N-1+itemp,i)=1;
                end
                if (B(ii1 ,jj1 ,i)==1)
                    AA(N-1+itemp,i)=1;
                end
                itemp=itemp+1;
            end
        end
    end
end

end
