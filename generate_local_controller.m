function K = generate_local_controller(edges,cost_list,alfa)
%Generate a local controller. Where each link finds the flow that would be
%optimal for the graph containing the source, and all the decendants of the
%source.
%Only used as initial guess for the "optimized" local controller.

%We assume that the graph is rooted. And thus each node has maximum one
%parent.
    nbr_of_nodes = size(edges,1)+1;
    nbr_of_edges = size(edges,1);
    W = zeros(nbr_of_nodes);
    inc_link_vec = zeros(1,nbr_of_nodes); %contains the number of the link that is incoming to a node (zero for the root).
    for i = 1:nbr_of_edges
        W(edges(i,1),edges(i,2)) = 1; %Used to find easily find all children for a node
        inc_link_vec(edges(i,2)) = i;
    end
    K = zeros(nbr_of_edges,nbr_of_nodes+nbr_of_edges);
    nbr_of_producers = 0; %Might change later.
    input_state_shift = nbr_of_producers + nbr_of_nodes;
    for i =1:nbr_of_edges %loop over all edges
        source = edges(i,1); %The source of the edge
        dest = edges(i,2); %and the desination
        other_children = setdiff(find(W(source,[1:end])),dest); %the other children of the source
        gamma_d = cost_list(dest)*alfa^2;
        gamma_u = 1./(sum(1./(alfa^2.*cost_list(other_children)))+1/cost_list(source));
        %Downstream set in the subgraph contains only the destination
        K(i,dest) = -alfa*gamma_d/(gamma_u+gamma_d);
        K(i,input_state_shift+i) = -alfa*gamma_d/(gamma_u+gamma_d); %we know parentlink has number i
        %Upstream set in the subgraph contains the source and all other
        %nodes:
        for j = other_children
            K(i,j) = alfa*gamma_u./(gamma_u+gamma_d);
            K(i,input_state_shift + inc_link_vec(j)) = alfa*gamma_u/(gamma_u+gamma_d);
        end
        K(i,source) = alfa*gamma_u/(gamma_u+gamma_d);
        if inc_link_vec(source)>0 %source might be root.
           K(i,input_state_shift+inc_link_vec(source) ) = alfa*gamma_u/(gamma_u+gamma_d);
        end
    end
end
        
       

