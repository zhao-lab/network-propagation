function x=opt_com_weight(H,i,N_p,com_list)    %take in a matrix as L^T*L, where L is the deviation vector 
%return a optimized communication weight, return the number of particles
%that i-th filter should collect
N=size(H,1);
lb=zeros(N,1);
lb(i)=0.25;    %such that the i-th filter at least has 10% remaining
ub=ones(N,1);
Aeq=ones(1,N);
beq=1;
not_com_list=setdiff(1:N,com_list);
for k=1:length(not_com_list)
    index=not_com_list(k);
    if index~=i
        Aeq_plus=zeros(1,N);
        Aeq_plus(index)=1;
        Aeq=[Aeq;Aeq_plus];
        beq=[beq;0];
    end
end
f=zeros(N,1);
y=quadprog(H,f,[],[],Aeq,beq,lb,ub);
x=round(y/y(i)*N_p); 
end

