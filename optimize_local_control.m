function [K_local_opt] = optimize_local_control(K_guess,edges,q_vec,alfa,sparsity_pattern)
%Searches for the locally optimal controller with presecribed sparsity
%only allows for q_i = 1
nbr_nodes = length(q_vec);
nbr_edges = length(edges);
nbr_states = nbr_nodes + nbr_edges;
[ A,B,Q,R ] = generate_graph(edges,[], q_vec,[] );
A = alfa*A;
if q_vec ~= ones(size(q_vec))
    error('Only supports q_vec = 1')
end
C = zeros(nbr_nodes,nbr_states);
C(1:nbr_nodes,1:nbr_nodes) = eye(nbr_nodes); %Only penalize Node states.

%Calculate the norm given a vectorized gain matrix
closed_loop_norm_vec = @(k) norm(ss(A+B*reformat(k),eye(size(A)),C,0,-1));
K_guess_vec = K_guess(:);
options = optimoptions('fmincon','MaxFunctionEvaluation',10000);
K_local_opt = reformat(fmincon(closed_loop_norm_vec,K_guess_vec(find(K_guess_vec)),[],[],[],[],[],[],[],options));


%Takes a vectorized version of the gain matrix K and convert it to a matrix
%of correct size
function [K] = reformat(k)
    entries = find(sparsity_pattern);
    K_vec = zeros(size(sparsity_pattern(:)));
    K_vec(entries) = k;
    K = reshape(K_vec,size(sparsity_pattern));
end

end




