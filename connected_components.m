function components = connected_components(A,N)
% A is the adjacency matrix
% N is the number of agents composing the connected components

n = size(A, 1);
nodes = 1:n;
components = {};
count = 0;

% Generate all combinations of N nodes
combs = nchoosek(nodes, N);

for i = 1:size(combs, 1)
    subset = combs(i, :);
    subgraph = A(subset,subset);
    G = graph(subgraph);
    bins = conncomp(G);
    if (max(bins)==1)
        count = count + 1;
        components{count} = subset;
    end
end

end



% % Example adjacency matrix from the user
% adj_matrix = [0 1 0 0 0 0 0 0 0 1;
%               1 0 1 1 0 0 0 0 0 0;
%               0 1 0 0 0 1 0 0 1 0;
%               0 1 0 0 0 1 1 0 1 0;
%               0 0 0 0 0 0 0 1 0 1;
%               0 0 1 1 0 0 0 0 0 0;
%               0 0 0 1 0 0 0 1 0 1;
%               0 0 0 0 1 0 1 0 1 0;
%               0 0 1 1 0 0 0 1 0 0;
%               1 0 0 0 1 0 1 0 0 0];
          
