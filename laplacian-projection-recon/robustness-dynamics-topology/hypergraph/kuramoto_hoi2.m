function f=kuramoto_hoi2(t,phi,EdgeList,TriangleList,omega,alpha,I)

[N,n1]=size(phi);

coup_rete=zeros(N,1);
coup_simplicial=zeros(N,1);

%pairwise terms
for ii=1:length(EdgeList)
    i1=EdgeList(ii,1);
    i2=EdgeList(ii,2);
    coup_rete(i1)=coup_rete(i1)+sin(phi(i2) - phi(i1));
    coup_rete(i2)=coup_rete(i2)+sin(phi(i1) - phi(i2));
end

%higher-order interactions
[mtrianglelist,ntrianglelist]=size(TriangleList);
for ii=1:mtrianglelist
    i1=TriangleList(ii,1);
    i2=TriangleList(ii,2);
    i3=TriangleList(ii,3);
    coup_simplicial(i1)=coup_simplicial(i1)+sin(phi(i2)+phi(i3) - 2 * phi(i1));
    coup_simplicial(i2)=coup_simplicial(i2)+sin(phi(i1)+phi(i3) - 2 * phi(i2));
    coup_simplicial(i3)=coup_simplicial(i3)+sin(phi(i1)+phi(i2) - 2 * phi(i3));
end


% Compute oscillator phase derivative.
f = omega + alpha(1)*coup_rete + alpha(2)*coup_simplicial+I;
end
