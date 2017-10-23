function initialize_pf_CMM(Np,Nsv,N,common_error,distance_dyn)
global pf;
for k=1:Np
    pf(k).common=common_error+0.1*randn(Nsv,1);  %common error stored as row vector
%     if distance_dyn > 2000
%         pf(k).weight=0.3/Np;
%     else
        pf(k).weight=1/Np;
%     end
end
for k=1:Np
%     pf(k).mu{1}=[16;v+0.3*rand(1);248.25;0];
%     pf(k).mu{2}=[473;-v+0.3*rand(1);251.75;0];
%     pf(k).mu{3}=[251.75;0;16;v+0.3*rand(1)];
%     pf(k).mu{4}=[248.25;0;473;-v+0.3*rand(1)];
    pf(k).mu{1}=[0;0;0;0];
    pf(k).mu{2}=[0;0;0;0];
    pf(k).mu{3}=[0;0;0;0];
    pf(k).mu{4}=[0;0;0;0];
    for j=5:N
        index=mod(j,4);
        if index==0
            index=4;
        end
        pf(k).mu{j}=pf(k).mu{index};
    end
end
for k=1:Np
    for j=1:N
    pf(k).cov{j}=eye(4)*10;
    end
end
end
