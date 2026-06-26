function AA = Tri_convert_to_AA(A,TriangleList)
N = length(A);              %网络中振子的数量
%number of unknowns for each node
H=(N-1)+(N-1)*(N-2)/2;
AA=zeros(H,N); %"true" matrix

for i=1:N
    AA(1:N-1,i)=A(i,[1:i-1 i+1:N]);
    itemp=1;
    for ii1=[1:i-1 i+1:N]
        for jj1=[1:i-1 i+1:N]
            if jj1>ii1
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
end

end