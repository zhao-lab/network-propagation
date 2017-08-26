function initialize_pf_CMM(Np,Nsv,N,common_error,v);
global pf;
for i=1:N
for k=1:Np
    pf{i}(k).common=common_error+0.1*randn(Nsv,1);  %common error stored as row vector
    pf{i}(k).weight=1/Np;
end
for k=1:Np
%     pf{i}(k).mu{1}=[16;v+0.3*rand(1);248.25;0];
%     pf{i}(k).mu{2}=[473;-v+0.3*rand(1);251.75;0];
%     pf{i}(k).mu{3}=[251.75;0;16;v+0.3*rand(1)];
%     pf{i}(k).mu{4}=[248.25;0;473;-v+0.3*rand(1)];
pf{i}(k).mu{1}=[0;0;0;0];
pf{i}(k).mu{2}=[0;0;0;0];
pf{i}(k).mu{3}=[0;0;0;0];
pf{i}(k).mu{4}=[0;0;0;0];
    for j=5:N
        index=mod(j,4);
        if index==0
            index=4;
        end
        pf{i}(k).mu{j}=pf{i}(k).mu{index};
    end
end
for k=1:Np
    for j=1:N
    pf{i}(k).cov{j}=eye(4)*10;
    end
end
end
end
