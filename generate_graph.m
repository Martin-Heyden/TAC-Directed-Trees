 function [ A,B,Q,R ] = generate_graph(edges,prod, Qvec, Rvec )
%Generate graph defined by edges. First node states, then delay
%states (Prod, internal) Order correspond to order in prod, then edges.
%Returns the state transition matrices A,B corresponding to edges and prod
%and the cost matrices Q,R corresponding to Q_vec and R_vec
number_of_prod = length(prod);
number_of_int = size(edges,1);
number_of_inputs = number_of_int+number_of_prod;
number_of_nodes = size(edges,1)+1;
number_of_states = number_of_nodes+ number_of_inputs;
A = zeros(number_of_states,number_of_states);
B = zeros(number_of_states,number_of_inputs);
Q = zeros(number_of_states,number_of_states);
R = zeros(number_of_inputs,number_of_inputs);

A(1:number_of_nodes,1:number_of_nodes) = eye(number_of_nodes,number_of_nodes); % node update dependent on node level
Q(1:number_of_nodes,1:number_of_nodes) = diag(Qvec); %node level penalties
R(1:number_of_prod,1:number_of_prod) = diag(Rvec);%production penalties
i = 1;
while i<=number_of_prod %loop over producers
   dest = prod(i);
   %Update matrices
   B(i+number_of_nodes,i) = 1; 
   A(dest,number_of_nodes+i) = 1;
   i = i +1;
end
j = 1;
while j<=number_of_int %loop over internal flows
    i = j+number_of_prod;
    source = edges(j,1);
    dest = edges(j,2);
    %Update matrices
    B(source,i) = -1;
    B(i+number_of_nodes,i) = 1;
    A(dest,number_of_nodes+i) = 1;
    j = j+1;
end



end