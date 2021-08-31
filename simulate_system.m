function [cost,x_mat] = simulate_system(edges,q_vec,K,alfa,T,x0,w_mat,save_states,prod,r_vec)
%Simulate system defined by edges, alfa  cost q_vec, prod and r_vec, with given
%controller K. Simulation runs for T units, or until cost is zero if T<0.
%The noise (if any) is supplied so it can be the same for different
%simulations. States are saved and returned if save_states>0. Cost is
%always calculated

%assume no production if not specified
if nargin == 8
    prod = [];
    r_vec = [];
end
%T = -1 means run until cost is zero
if T<0 &  save_states>0
    error('SimulateSystem:invalidSettings','Can not save states (save_states>0) when horizon is not specified (T<0)');
end

%Define state space matrice
[ A,B,Q,R ] = generate_graph(edges,prod, q_vec,r_vec );
A = alfa*A;
nbr_states = length(A);

if save_states %structure for saving states
    x_mat = zeros(nbr_states,T+1);
    x_mat(:,1) = x0;
else
    x_mat = [];
end


go_on = 1;
t = 0;
x = x0;
agg_cost = x0'*Q*x0;

while go_on
   t = t+1;
   u = K*x;%calculate input
   %update states
   if w_mat == -1 %no noise
       x = A*x+B*u;
   else %with noise
    x = A*x + B*u + w_mat(:,t);
   end
   
   if save_states %save states, if they are to be saved
       x_mat(:,t+1) = x;
   end
   if length(r_vec)==0 %if no input penalty
        cost = x'*Q*x;
   else %otherwise there is a input penalty
       cost = x'*Q*x + u'*R*u;
   end
   agg_cost = agg_cost + cost; %update cost
   if T<0 %If stopping based on cost
       if cost <1e-6
           go_on = false;
       end
   else %If stopping based on time
       if t == T
           go_on=false;
       end
   end   
end

end

