function [mu,cov]=pf_to_gaussian_RBPF(j)
global pf;
mu=0;
%N=length(pf{j});
N = 50;
for k=1:N
    mu=mu+pf{j}(k).weight*pf{j}(k).mu{j}([1,3]);
end
cov=zeros(2,2);
for k=1:N
    cov=cov+pf{j}(k).weight*(pf{j}(k).cov{j}([1,3],[1,3])+(pf{j}(k).mu{j}([1,3])-mu)*(pf{j}(k).mu{j}([1,3])-mu)');
end
end
