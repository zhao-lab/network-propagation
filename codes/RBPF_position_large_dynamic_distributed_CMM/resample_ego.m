function resample_ego
global pf_pose;
pf_copy=pf_pose;
% do some stuff here
size_samples=pf_pose.N;
pf_pose.weight=pf_pose.weight/sum(pf_pose.weight);

c(1)=pf_pose.weight(1);
for i=2:size_samples
    c(i)=c(i-1)+pf_pose.weight(i);
end
i=1;
u(1)=1/size_samples*rand(1);
for j=1:size_samples
    while u(j) > c(i)
        i=i+1;
    end
    pf_pose.state(j,:)=pf_copy.state(i,:);
    pf_pose.mp(j,:)=pf_copy.mp(i,:);
    pf_pose.mpvar(j,:)=pf_copy.mpvar(i,:);
    u(j+1)=u(j)+1/size_samples;
end
pf_pose.weight=1/pf_pose.N*ones(pf_pose.N,1);
end
