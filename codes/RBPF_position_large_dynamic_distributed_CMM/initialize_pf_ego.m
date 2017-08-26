function initialize_pf_ego(Np_pose,state,derivative,sv_number)  %state contains x,y and clock bias 
%equivalent distance. sv_number gives the number of satellite, other input
%argument should be added later. e.g., satellite id...
global pf_pose;
for k=1:Np_pose
    pf_pose.state(k,:)=[state, derivative];   %state contains 6 scalar, for higher efficiency, consider putting derivative into Kalman filter
end
pf_pose.mp=0.1*randn(Np_pose,sv_number);
pf_pose.mpvar=ones(Np_pose,sv_number);
pf_pose.weight=1/Np_pose*ones(Np_pose,1);
pf_pose.N=Np_pose;
pf_pose.weight_time=[];
end

    