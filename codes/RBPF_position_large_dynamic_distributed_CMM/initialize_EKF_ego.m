function initialize_EKF_ego(state,derivative,sv_number)  %state contains x,y and clock bias 
%equivalent distance. sv_number gives the number of satellite, other input
%argument should be added later. e.g., satellite id...
global EKF;
EKF.N=6+2*sv_number;   %length of state vecotr
EKF.mu=[state(1);derivative(1);state(2);derivative(2);state(3);derivative(3)];
EKF.mu=[EKF.mu;zeros(2*sv_number,1)];  %initialize multipath with 0
EKF.cov=eye(EKF.N,EKF.N);   %initialize covariance
EKF.cov(7:end,7:end)=100*EKF.cov(7:end,7:end);  %increase the covariance of the multipath error
end