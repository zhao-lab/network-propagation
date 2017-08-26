function pf_update_ego(pr,svxyzmat,sigma_thermal2,orgxyz);
global pf_pose;
sv_n=length(pr);
for k=1:pf_pose.N
    pose_local=[pf_pose.state(k,1:2),0];  %here, the z coordinate is given 0
    pose_xyz=enu2xyz(pose_local,orgxyz);  %give xyz as column vector
    for j=1:sv_n
        if j<=2
    range_pre(k,j)=norm(pose_xyz.'-svxyzmat(j,:))+pf_pose.state(k,3)+pf_pose.mp(k,j);  %prediction of pr measurement, distance+clock bias+multipath
        else
    range_pre(k,j)=norm(pose_xyz.'-svxyzmat(j,:))+pf_pose.state(k,3);
        end
        
    innovation(k,j)=pr(j)-range_pre(k,j);
    
    if j<=2
    Kalman_gain(k,j)=pf_pose.mpvar(k,j)/(pf_pose.mpvar(k,j)+sigma_thermal2);
    end
    
    end
end
% It should be noted that only multipath gets updated, while the pose is
% not. This is quite inefficient. Better scheme can be used which is
% similar to FASTslam2.0
pf_pose.mp=pf_pose.mp+Kalman_gain.*innovation(:,1:2);
pf_pose.mpvar=(1-Kalman_gain).*pf_pose.mpvar;   %Kalman update for multipath
% variance of multipath decreases monotonically, this should be addressed.
% Perhaps by including the derivative of multipath error

%calculate particle weight
Q=[pf_pose.mpvar+sigma_thermal2,ones(pf_pose.N,sv_n-2)*sigma_thermal2];

%Q=zeros(500,4)+sigma_thermal2;
for k=1:pf_pose.N
    for j=1:sv_n
     w=1/sqrt(2*pi*Q(k,j));
     w=w*exp(-0.5*innovation(k,j)^2/Q(k,j));
     pf_pose.weight(k)=pf_pose.weight(k)*w;
    end
end
pf_pose.weight_time=[pf_pose.weight_time,pf_pose.weight];
end


    