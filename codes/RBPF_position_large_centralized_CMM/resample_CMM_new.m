function resample_CMM_new(m,distance_dyn,confi_rang,num_confi,num_block_drop)
global pf;
pf_copy=pf{m};
% do some stuff here
% size_samples = Np in main function.
size_samples=length(pf{m});
S=0;
for k=1:size_samples
    S=S+pf{m}(k).weight;
end
for k=1:size_samples
    pf{m}(k).weight=pf{m}(k).weight/S;
end

c(1)=pf{m}(1).weight;
for i=2:size_samples
    c(i)=c(i-1)+pf{m}(i).weight;
end
i=1;
u(1)=1/size_samples*rand(1);
for j=1:size_samples
    while u(j) > c(i)
        i=i+1;
    end
    pf{m}(j)=pf_copy(i);
    u(j+1)=u(j)+1/size_samples;
end

times = 1;
% ratio_in_confi = num_conf/size_sample;
rate_of_sample = 1/(size_samples + (times-1)*num_confi - num_block_drop);

for k=1:size_samples
%     if equal: pf(k).weight=1/size_samples;
    if distance_dyn < confi_rang
        pf{m}(k).weight=times*rate_of_sample;  %*pf(k).weight;
    else
        pf{m}(k).weight=rate_of_sample;  %*pf(k).weight;
    end
end

% for k=1:size_samples
%     pf(k).weight=1/size_samples;
%     if distance_dyn{k} > confi_rang
%         pf(k).weight=0.2*pf(k).weight;
%     end
% end

% for k=1:size_samples
%     pf(k).weight=1/size_samples;
%     if distance_dyn <=20
%         pf(k).weight = 0.2*pf(k).weight;
%     elseif distance_dyn > confi_rang
%         pf(k).weight=0.2*pf(k).weight;
%     end
% end


end
