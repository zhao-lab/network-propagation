function pf_prediction_ego(sigma_x,sigma_x_dot,sigma_b,sigma_b_dot,deltat)
global pf_pose;

pf_pose.state(:,1)=pf_pose.state(:,1)+pf_pose.state(:,4)*deltat;
pf_pose.state(:,2)=pf_pose.state(:,2)+pf_pose.state(:,5)*deltat;
pf_pose.state(:,3)=pf_pose.state(:,3)+pf_pose.state(:,6)*deltat;

s_y=0.01;
crossx=sqrt(sigma_x*sigma_x_dot);
crossb=sqrt(sigma_b*sigma_b_dot);
cov=[sigma_x,  0,  0, crossx,  0,  0;
     0,  s_y*sigma_x,  0,  0,    s_y*crossx, 0;
     0,     0,    sigma_b, 0,  0,   crossb;
     crossx, 0,    0,   sigma_x_dot, 0, 0;
     0,  s_y*crossx,   0,       0, s_y*sigma_x_dot, 0;
     0,    0,     crossb,  0,   0, sigma_b_dot];
sample=mvnrnd([0,0,0,0,0,0],cov,pf_pose.N);
pf_pose.state=pf_pose.state+sample;
pf_pose.mpvar=pf_pose.mpvar+0.25; %add variance to model change in mpvar
end
