function [edges] = generate_edge_list(depth,nbr_of_children)
%Generate the edge list for a graph with depth layers, and where each node
%(except those in the last layer) has nbr_children children.
% The root is node 1. the nodes with depth 2 are 2... 1+nbr_of_children
%and so on

edges = [];
next_node = 1;
next_child = 2;
for i = 1:depth-1
    for j = 1:nbr_of_children^(i-1)%nbr of nodes in layer.
        for k = 1:nbr_of_children %nbr of children for each node in layer
            edges = [edges; [next_node,next_child]]; 
            next_child = next_child+1;
        end
        next_node = next_node + 1;
    end         
end
end

