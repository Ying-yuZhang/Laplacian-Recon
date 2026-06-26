function A=erdos_reyni(n,p)
% ERDOS_REYNI Generates a random Erdos-Reyni (Gnp) graph
%
% A=erdos_reyni(n,p) generates a random Gnp graph with n vertices and where
% the probability of each edge is p.  The resulting graph is symmetric.
%

% Example:
%   A = erdos_reyni(100,0.05);

   A = triu(rand(n)<p,1);
    A = sparse(A);
    A = A+A';
    
    
% C=0;
% while(max(C)~=n)
%     clear A
%     A = triu(rand(n)<p,1);
%     A = sparse(A);
%     A = A+A';
%     [aa,C] = components(A);
% end
