function weight_CMM(vehicle_id,pf_id,ang)
global pf;
N_mc=30;  % #of sample for Monte Carlo integration
Np=length(pf{pf_id});
road_normal=[sin(ang);-cos(ang)];
for k=1:Np
    sample(:,1:2)=mvnrnd([pf{pf_id}(k).mu{vehicle_id}(1),pf{pf_id}(k).mu{vehicle_id}(3)],pf{pf_id}(k).cov{vehicle_id}([1,3],[1,3]),N_mc);
% x=(sample(:,1)>247).*(sample(:,1)<253);
% y=(sample(:,2)>247).*(sample(:,2)<253);
% pf{pf_id}(k).weight=pf{pf_id}(k).weight*sum(x|y);

pf{pf_id}(k).weight=pf{pf_id}(k).weight*sum((sample*road_normal<1));  %assume that this project larger than 1 the road constraint would be active
end
end
