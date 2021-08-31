%Code for verifying theorem 1 and 2 for rooted and non rooted trees
%% Test rooted Tree (edges can be changed to represent any rooted tree)
edges = generate_edge_list(5,2);
nbr_edges = size(edges,1);
qvec = rand(length(edges)+1,1);
decay = 0.95;
r = 10

%Calculate controller using theorem 1 and 2
K = synthesis_rooted_tree(edges,qvec,decay,r); 


[ A,B,Q,R] = generate_graph(edges, [1], qvec,[r] ); %One producer
A = decay*A;
[X,L,G] = dare(A,B,Q,R); %Calculate controller numerically
K_opt = -G;

norm(K-K_opt,Inf) %Should be ~zero

%% Test non rooted tree (keep edges as is, or update downstream and upstream manualy) 

edges = [4 1;5 1; 5 2;6 3;7 5;7 6;8 6; 10 7;11 8;11 9]; %corresponds to Fig 6 + one additional node
nbr_nodes = length(edges)+1;
nbr_edges = length(edges);
qvec = rand(nbr_nodes,1);
upstream = cell(1,nbr_edges); %upstream set for each each
downstream = cell(1,nbr_edges); %Downstream set for each each

alfa = 0.95;

%%% This part is manual %%%
upstream{1} = [4]; downstream{1} = [1 2 5];
upstream{2} = [2 5]; downstream{2} = [1 4];
upstream{3} = [1 4 5]; downstream{3} = [2];
upstream{4} = [6]; downstream{4} = [3];
upstream{5} = [3 6 7 8]; downstream{5} = [1 2 4 5];
upstream{6} = [1 2 4 5 7]; downstream{6} = [3 6 8];
upstream{7} = [8]; downstream{7} = [1 2 3 4 5 6 7];
upstream{8} = [10]; downstream{8} = [1 2 3 4 5 6 7 8 9 11];
upstream{9} = [9 11]; downstream{9} = [1 2 3 4 5 6 7 8 10];
upstream{10} = [1 2 3 4 5 6 7 8 10 11]; downstream{10} = [9];
parent_link = { [1 2] 3 4 [] 5 [6 7] 8 9 10 [] [] }
depth = [3 3 3 2 2 2 1 1 1 0 0]; %depth of each node
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

K = zeros(nbr_edges,nbr_nodes); 
%Optimal internal flows:
for i = 1:nbr_edges
    source_depth = depth(edges(i,1));
    gamma_u = 1/(sum(1./(alfa.^(2*depth(upstream{i})-2*source_depth).*qvec(upstream{i})')));
    gamma_d = 1/(sum(1./(alfa.^(2*depth(downstream{i})-2*source_depth).*qvec(downstream{i})')));
    K(i,[upstream{i} nbr_nodes+[parent_link{upstream{i}}]]) = alfa*gamma_u/(gamma_u+gamma_d);
    K(i,[downstream{i} nbr_nodes+[parent_link{downstream{i}}]]) = - alfa* gamma_d/(gamma_u+gamma_d);
end
if 0 %no production
    prod = []
    r_vec = []
else %With production
    prod = [10 11]
    r_vec = rand(1,2)
    %Optimal production:
    gamma_v = 1/(sum(1./(alfa.^(2*depth).*qvec')));
    R = 1/(sum(1./r_vec));
    X = -1/2*((1-alfa^2)*R-alfa^2*gamma_v) + sqrt(alfa^2*gamma_v*R+1/4*((1-alfa^2)*R-alfa^2*gamma_v)^2);
    K = [-ones(1,nbr_edges*2+3)*R/r_vec(1)*X/(X+R)*alfa;
        -ones(1,nbr_edges*2+3)*R/r_vec(2)*X/(X+R)*alfa;
        [K(:,1:nbr_nodes), zeros(nbr_edges,2), K(:,nbr_nodes+1:end)]
        ];
    %adding dependence on delayed production for internal flows
    for i = [8 9 10]+2
            K(i,nbr_nodes+1) = K(i,10);
            K(i,nbr_nodes+2) = K(i,11);
    end
end
   
%Compare to numerical solution
[ A,B,Q,R] = generate_graph(edges, prod, qvec,r_vec ); %One producer
A = alfa*A;
[X,L,G] = dare(A,B,Q,R);
K_opt = -G;

norm(K-K_opt,Inf)

