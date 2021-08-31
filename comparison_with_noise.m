%Code for generating Fig 11.
rng(123)
%%% Settings %%%
if 0 % String graph
    max_depth = 100;
    start_depth = 5;
    depth_step = 5;
    nbr_children = 1;
else %Binary tree
    max_depth = 8;
    start_depth = 2;
    depth_step = 1;
    nbr_children = 2;
end
depth_vec = start_depth:depth_step:max_depth;
nbr_depth = length(depth_vec)
nbr_simulation = 1000; 
decay = 0.99;
cost_threshold = 1e-6;
T = 100

%%% Data Matrices %%%
opt_cost = zeros(1,nbr_depth);
local_cost = zeros(1,nbr_depth);
opt_prod_cost = zeros(1,nbr_depth);

for ind  = 1:nbr_depth
    depth = depth_vec(ind)
    %%% Generate System for current depth %%%
    edges = generate_edge_list(depth,nbr_children);
    nbr_nodes = size(edges,1)+1;
    nbr_edges = size(edges,1);
    qvec = ones(nbr_nodes,1);
    %Generate the different controllers
    K_opt = synthesis_rooted_tree(edges,qvec,decay);
    K_opt_prod = synthesis_rooted_tree(edges,qvec,decay,10);
    K_local_init = generate_local_controller(edges,qvec,decay);
    K_local = optimize_local_control(K_local_init,edges,qvec,decay,K_local_init);
    for i = 1:nbr_simulation
        x0 = zeros(length(K_opt),1); %Same initial condions for all cases
        w_mat = randn(length(K_opt),T)*0.01; %And same noise for all cases
        %Simulate the different systems
        [cost_opt,~] = simulate_system(edges,qvec,K_opt,decay,T,x0,w_mat,0);
        [cost_opt_prod,~] = simulate_system(edges,qvec,K_opt_prod,decay,T,[x0;0],[zeros(1,T);w_mat],0,1,0.1);
        [cost_loc,~] = simulate_system(edges,qvec,K_local,decay,T,x0,w_mat,0);
        %Add the cost
        opt_cost(ind) = opt_cost(ind) + cost_opt;
        local_cost(ind) = local_cost(ind) + cost_loc;
        opt_prod_cost(ind) = opt_prod_cost(ind) + cost_opt_prod;
    end
    %Normalize by number of nodes
    opt_cost(ind) = opt_cost(ind)/nbr_nodes;
    opt_prod_cost(ind) = opt_prod_cost(ind)/nbr_nodes;
    local_cost(ind) = local_cost(ind)./nbr_nodes;
end



%%% Plotting %%%
hold off
plot(depth_vec,local_cost,'-x','linewidth',2);
hold on
plot(depth_vec,opt_cost,'-x','linewidth',2);
plot(depth_vec,opt_prod_cost,'-d','linewidth',2);
legend('Local Controller','Optimal Controller','Optimal with Production','Location','NorthEast')
xlabel('Network Depth','fontsize',12)
ylabel('Cost Per Node','fontsize',12)
set(gcf,'position',[900,500,550,250]) %x0 y0, width height
axis([start_depth-1,max_depth+1,0,max(local_cost)*1.1])
