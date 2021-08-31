function K = synthesis_rooted_tree(edges,qvec,decay,r)
% Find the optimal controller for the system defined by edges and qvec.
% Assumes that edges is a rooted tree and that when a node is defined (by
% appearing as a destination in edges) then all its children will be bellow
% it in the list.
% r negative means no production. Assumes first node is producer (if any)

%%% The implementation is not the most efficient (especially not in matlab)
%%% However, it emphezise how the synthesis can be made in a local way

if nargin == 3
    r = -1 %Corresponds to no input (there can only be one for rooted trees)
end

%Data structures
nbr_edges = size(edges,1);
nbr_nodes = nbr_edges + 1;

K = zeros(nbr_edges,nbr_nodes+nbr_edges);
incoming_link_list = zeros(nbr_nodes,1); % A list of the unique incoming link for every node
children_sets = cell(nbr_nodes,1);
decendent_sets = cell(nbr_nodes,1); %the decendents for every node
aggregate_cost = qvec; %initializes the nodes with no children correctly.

nbr_of_nodes = size(edges,1)+1;

%first pass, loop over edges, starting at "the bottom"
for i = nbr_edges:-1:1
    source = edges(i,1);
    dest = edges(i,2);
    incoming_link_list(dest) = i;
    %Children set is old children set + the new child
    children_sets{source} = [children_sets{source},dest];
    %Decendent set is the old decendent set, the new child, and the
    %decendents of the new child.
    decendent_sets{source} =  [decendent_sets{source}, dest,decendent_sets{dest}];
    %Update aggregate cost of dest as in paper (all its children has been processed now).
    aggregate_cost(dest) = inv(1/qvec(dest) + sum(1./(decay^2*aggregate_cost(children_sets{dest}))));
    
    %We are now ready to find the optimal flow for every outgoing link from
    %the destination, as all children and grandchildren has been processed
    parent = dest; %Makes this section easier to read.
    for child = children_sets{parent}%loop over children          
        link_number = incoming_link_list(child); %Which gives all outgoing links
        downstream_set = [child,decendent_sets{child}]; %find the downstream
        upstream_set = setdiff([parent,decendent_sets{parent}],downstream_set);%and upstream set
        gamma_d = decay^2*aggregate_cost(child);%Find gamma_d
        gamma_u = inv(1/aggregate_cost(parent)-1/(decay^2*aggregate_cost(child)));%and then calculate gamma_u
        K(link_number,downstream_set) = -decay*gamma_d/(gamma_u+gamma_d);%dependence of node levels in downstream set
        K(link_number,nbr_of_nodes + incoming_link_list(downstream_set)) = -decay*gamma_d/(gamma_u+gamma_d);%transportation states in downstream set
        K(link_number,upstream_set) = decay*gamma_u/(gamma_u+gamma_d);%node levels in upstream set
        K(link_number,nbr_of_nodes + incoming_link_list(upstream_set)) =decay*gamma_u/(gamma_u+gamma_d);%transportation states in upstream set
    end
end
%The first node will not have a parent, and is thus not covered by the foor
%loop above.
aggregate_cost(1) = inv(1/qvec(1) + sum(1./(decay^2*aggregate_cost(children_sets{1}))));
parent = 1; 
for child = children_sets{parent}           
    link_number = incoming_link_list(child);
    downstream_set = [child,decendent_sets{child}];
    upstream_set = setdiff([parent,decendent_sets{parent}],downstream_set);
    gamma_d = decay^2*aggregate_cost(child);
    gamma_u = inv(1/aggregate_cost(parent)-1/(decay^2*aggregate_cost(child)));
    K(link_number,downstream_set) = -decay*gamma_d/(gamma_u+gamma_d);
    K(link_number,nbr_of_nodes + incoming_link_list(downstream_set)) = -decay*gamma_d/(gamma_u+gamma_d);
    K(link_number,upstream_set) = decay*gamma_u/(gamma_u+gamma_d);
    %NOTE: incoming_link to node 1 is zero, indicate that there is no link
    K(link_number,nbr_of_nodes + incoming_link_list(upstream_set(2:end))) =decay*gamma_u/(gamma_u+gamma_d);
end
if r>0 %if there is production
    gamma_v = aggregate_cost(1); %aggregate cost for entire graph
    %Corresponds to a discrete time Riccati equation:
    X = -1/2*((1-decay^2)*r-decay^2*gamma_v) + sqrt(decay^2*gamma_v*r+1/4*((1-decay^2)*r-decay^2*gamma_v)^2);
    %Add a new row to K
    K = [ones(1,nbr_edges*2+2)*-X/(X+r)*decay;
        [K(:,1:nbr_nodes), zeros(nbr_edges,1), K(:,nbr_nodes+1:end)]
        ];
    %Update the old rows that are affected by the production (links going
    %out of node 1)
    for i = decendent_sets{1}
        K(incoming_link_list(i)+1,nbr_nodes+1) = K(i,1); 
    end
end




end

