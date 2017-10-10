function weight_CMM(vehicle_id,ang,distance)
global pf;
N_mc=100;  % #of sample for Monte Carlo integration
Np=length(pf);
road_normal=[sin(ang);-cos(ang)];
for k=1:Np
    sample(:,1:2)=mvnrnd([pf(k).mu{vehicle_id}(1), ...
        pf(k).mu{vehicle_id}(3)],pf(k).cov{vehicle_id}([1,3],[1,3]),N_mc);
    % x=(sample(:,1)>247).*(sample(:,1)<253);
    % y=(sample(:,2)>247).*(sample(:,2)<253);
    if distance > 3000
%        pf(k).weight=pf(k).weight*0.2*sum((sample*road_normal<1));
        pf(k).weight=0.017;
    else
        pf(k).weight=1.2*pf(k).weight*sum((sample*road_normal<1));
    end
    %pf(k).weight=pf(k).weight*sum(x|y);
end
end
