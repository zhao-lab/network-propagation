function y=pf_get_common(k,H)
global pf;
N=size(H,1);
mu=zeros(N,1);
for j=1:length(pf{k})
mu=mu+pf{k}(j).common;
end
mu=inv(H.'*H)*H.'*mu/N;
y=mu(1:2);
end