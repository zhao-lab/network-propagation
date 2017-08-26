function resample_CMM(pf_id)
global pf;
pf_copy=pf{pf_id};
% do some stuff here
size_samples=length(pf{pf_id});
S=0;
for k=1:size_samples
    S=S+pf{pf_id}(k).weight;
end
for k=1:size_samples
    pf{pf_id}(k).weight=pf{pf_id}(k).weight/S;
end

c(1)=pf{pf_id}(1).weight;
for i=2:size_samples
    c(i)=c(i-1)+pf{pf_id}(i).weight;
end
i=1;
u(1)=1/size_samples*rand(1);
for j=1:size_samples
    while u(j) > c(i)
        i=i+1;
    end
    pf{pf_id}(j)=pf_copy(i);
    u(j+1)=u(j)+1/size_samples;
end
for k=1:size_samples
    pf{pf_id}(k).weight=1/size_samples;
end
end