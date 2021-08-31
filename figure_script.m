%Used to generate Figure 10

rng(123)

edges = generate_edge_list(20,1);
plot_ind = [1:5:16,20]
T = 30

nbr_edges = length(edges);
nbr_nodes = length(edges)+1;
q_vec = ones(1,nbr_nodes);
alfa =0.99;

w_mat = -1; %no noise
save_states = 1;

%Generate controllers
K_loc = generate_local_controller(edges,q_vec,alfa);
K_loc_opt = optimize_local_control(K_loc,edges,q_vec,alfa,K_loc);
K_opt = synthesis_rooted_tree(edges,q_vec,alfa)


n = nbr_nodes*2-1
x0 = mvnrnd(zeros(n,1),(n/(n-1))*(eye(n)-1/n*ones(n)))'; %Same initial conditions for both cases

%Carry out simulations
[cost_loc, x_mat_loc] = simulate_system(edges,q_vec,K_loc_opt,alfa,T,x0,w_mat,save_states);
[cost_opt, x_mat_opt] = simulate_system(edges,q_vec,K_opt,alfa,T,x0,w_mat,save_states);

%Plotting
figure(1)
clf
subplot(2,1,1)
hold on
for i =1:nbr_nodes
    if find(plot_ind == i),
        p = plot([0:1:T],x_mat_opt(i,1:T+1),'LineWidth',3) ;
    end
end
ax = gca;randn(nbr_nodes,1);
ax.FontSize = 10; 
box(ax,'on')
xlabel('Time (samples)','fontsize',12)
ylabel('$z_i[t]$','fontsize',12,'Interpreter','latex')
title('Optimal Controller', 'fontsize', 14)
legend(string(num2cell(abs(plot_ind-1))))
subplot(2,1,2)
hold on
for i =1:nbr_nodes
    if find(plot_ind == i),
        plot([0:1:T],x_mat_loc(i,1:T+1),'LineWidth',3) ;
    end
end
ax = gca;
ax.FontSize = 10; 
box(ax,'on')
xlabel('Time (samples)','fontsize',12)
ylabel('$z_i[t]$','fontsize',12,'Interpreter','latex')
title('Local Controller', 'fontsize', 14)
legend(string(num2cell(plot_ind-1)))
set(gcf,'position',[900,500,550,500]) %x0 y0, width height
