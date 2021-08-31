%Code for generating Figure 9
rng(456)
%%% Settings %%%
if 0 %string graph
    max_depth = 100;
    start_depth = 5;
    depth_step = 5;
    nbr_children = 1;
else %tree graph
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

%%% Data Matrices %%%4000
opt_cost = zeros(1,nbr_depth);
local_cost = zeros(1,nbr_depth);

for ind  = 1:nbr_depth
    depth = depth_vec(ind)
    %%% Generate System for current depth %%%
    edges = generate_edge_list(depth,nbr_children);
    nbr_nodes = size(edges,1)+1;
    nbr_edges = size(edges,1);
    qvec = ones(nbr_nodes,1);
    [ A,B,Q,~] = generate_graph(edges,[], qvec, [] ); %no production
    A = decay*A;
    %Generate controllers
    K_opt = synthesis_rooted_tree(edges,qvec,decay);
    K_local_init = generate_local_controller(edges,qvec,decay); %initial guess
    K_local = optimize_local_control(K_local_init,edges,qvec,decay,K_local_init);
    for i = 1:nbr_simulation
        %%% Generate Initial Conditions %%%
        n = nbr_nodes*2-1;
        x0 = mvnrnd(zeros(n,1),(n/(n-1))*(eye(n)-1/n*ones(n)))'; %Makes sum zero, and variance constant
        %%% Synthesize and simulate optimal controller %%%       
        opt_cost(ind) = opt_cost(ind) + simulate_controller(A,B,Q,K_opt,x0,cost_threshold)/nbr_simulation;
        local_cost(ind) = local_cost(ind) + simulate_controller(A,B,Q,K_local,x0,cost_threshold)/nbr_simulation; 
    end
    %Normalize by number of nodes
    opt_cost(ind) = opt_cost(ind)/nbr_nodes;
    local_cost(ind) = local_cost(ind)./nbr_nodes;
end



%%% Plotting %%%
hold off
plot(depth_vec,opt_cost,'-x','linewidth',2);
hold on
plot(depth_vec,local_cost,'-d','linewidth',2);
legend('Optimal Controller', 'Local Controller','Location','NorthWest')
xlabel('Network Depth','fontsize',12)
ylabel('Cost Per Node','fontsize',12)
set(gcf,'position',[900,500,550,250]) %x0 y0, width height

function total_cost = simulate_controller(A,B,Q,K,x0,cost_threshold) %Requires matlab 2016b or higher
    current_cost = intmax;
    x = x0;
    total_cost = 0;
    while current_cost>cost_threshold
        u = K*x;
        x = A*x + B*u;
        current_cost = x'*Q*x;
        total_cost = total_cost + current_cost;
    end
end