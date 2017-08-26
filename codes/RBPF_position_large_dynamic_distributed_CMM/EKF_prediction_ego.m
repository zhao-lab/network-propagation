function EKF_prediction_ego(sigma_a2,sigma_b2,sigma_d2,sigma_m2,sigma_md2,deltat)
%sigma_x and sigma_x_dot are the variance of unmodelled acceleration variance alone the lane
%sigma_b,sigma_b_dot are variance of clock bias
global EKF;
s_y=0.01;   %scaling factor for the lateral acceleration
Q_a=[sigma_a2*0.25*deltat^4, sigma_a2*0.5*deltat^3;
    sigma_a2*0.5*deltat^3,  sigma_a2*deltat^2];
Q_c=[sigma_b2*deltat^2+0.25*sigma_d2*deltat^4,0.5*sigma_d2*deltat^3;
    0.5*sigma_d2*deltat^3, sigma_d2*deltat^2];
Q_m=[sigma_m2*deltat^2+0.25*sigma_md2*deltat^4,0.5*sigma_md2*deltat^3;
    0.5*sigma_md2*deltat^3, sigma_md2*deltat^2];   %Noted that the sigma_a2 is of higher order in deltat, this is to account for the truncation error in the propagation equation
Sigma(1:2,1:2)=Q_a;  %lanewise
Sigma(3:4,3:4)=s_y*Q_a;  %lateral
Sigma(5:6,5:6)=Q_c;
A=[1 deltat;
    0 1];
for k=1:(EKF.N-6)/2
    Sigma(6+2*k-1:6+2*k,6+2*k-1:6+2*k)=Q_m;
end
%prediction
for k=1:EKF.N/2
    EKF.mu(2*k-1)=EKF.mu(2*k-1)+deltat*EKF.mu(2*k);
    EKF.cov(2*k-1:2*k,2*k-1:2*k)=A*EKF.cov(2*k-1:2*k,2*k-1:2*k)*A';
end
    EKF.cov=EKF.cov+Sigma; 
end