%this is a test program for particle filter mitigation of the multipath
%error on the state estimation. The methodology is to use particles to
%represent the vehicle state while every particle maintains a bank of
%Kalman filter for multipath. For every particle there are almost 24 multipath values, which 
% are assumed to be conditionally independent given vehicle pose. This
% assumption has to be checked later. Cooperative map matching should be
% added later as additional measurement to improve localization accuracy.

%created by macshen
clear all
close all
global pf_pose;
N=1; % # of total vehicles 
% load('distance_and_gradient.mat');  %contains grid_size, distance and gradient of distance from road centers
% global dis;
% global grad;   %global quantities to share with sub-function for interpolating distances
% dis=distance;    
% grad=grad_dis;

lane_width=3.5;   %width of a single lane
vehicle_width=1.8;  %use vehicle width to increase the positioning accuracy, it is assumed that the GPS receiver locates at the center of the vehicle
velocity=4;  %vehicle velocity
Ns=150;  % #of simulation time points
Np=500; % #of particles
Np_pose=100; % #of particles for ego state estimation and multipath mitigation
mpmat=cell(N);   %multipath error time history for N vehicles
for i=1:N
mpmat{i}=mpgen(24,1200,1,floor(10000*rand(1))+54321+i); %generate multipath error, use different random seeds
end

%mpmat{1}(1:Ns,1)=10000*sin(1/100*(1:Ns))';

orgllh = [40*pi/180 80*pi/180 0];  %origin of the local coordinate system, first latitude, second longitudinate 
orgxyz = llh2xyz(orgllh);   %convert the origin to the xyz ECEF cooridinate system
loadgps
startt=500; t = startt; deltat=0.1;
%segp = [150 90 .2; 150 90 .2; 150 90 .2];
%usrenu = pathgen([0 0 0],[5 0],segp,deltat);
usrenu{1}(1:Ns,1)=(6:velocity*deltat:6+velocity*deltat*(Ns-1))';
usrenu{1}(1:Ns,2)=50-lane_width/2;
usrenu{1}(1:Ns,3)=0;
% for k=1:3
% usrenu{k}(1:Ns,1)=(6:velocity:6+velocity*(Ns-1))';
% usrenu{k}(1:Ns,2)=50-lane_width/2;
% usrenu{k}(1:Ns,3)=0;
% end
% for k=4:6
% usrenu{k}(1:Ns,1)=(93:-velocity:93-velocity*(Ns-1))';
% usrenu{k}(1:Ns,2)=50+lane_width/2;
% usrenu{k}(1:Ns,3)=0;
% end
% for k=7:9
% usrenu{k}(1:Ns,2)=(6:velocity:6+velocity*(Ns-1))';
% usrenu{k}(1:Ns,1)=50+lane_width/2;
% usrenu{k}(1:Ns,3)=0;
% end
% for k=10:12
% usrenu{k}(1:Ns,2)=(93:-velocity:93-velocity*(Ns-1))';
% usrenu{k}(1:Ns,1)=50-lane_width/2;
% usrenu{k}(1:Ns,3)=0;             %generate the vehicles' path in the local coordinate system
% end
%EndLoop = max(size(usrenu));
%bar1 = waitbar(0,'Generating User Solutions...  ');
%id=[]; %added by macshen
pr_err=[]; %added by macshen
for i = 1:Ns  %draw the first step
    t = t + deltat;
    noise_scaling_factor=[0.5 0 0 0.5 0];
     
    
    for k=1:N   %GPS measurement loop for all the vehicles
    usrxyz{k}=enu2xyz(usrenu{k}(i,:),orgxyz);
    [svxyzmat{k},svid{k}] = gensv(usrxyz{k},t,0,orgxyz,7);   % here 5 is the mask angle to determine which satellites are visible to the vehicles, macshen
   % id=[id;svid(1)]; The last argument is to specify the wanted number of satellite and give the first group of  available satellites %added by macshen
    [prvec{k},adrvec{k},pr_err_vec{k}] = genrng(k,usrxyz{k},svxyzmat{k},svid{k},t,noise_scaling_factor,[],mpmat{k});  %pr_err_vec added by macshen, receiver identification k seems to be redundant here. the magnitude of noise flag is thermal,Tropospheric,SA,Multipath,Ionospheric
     pr_err=[pr_err;pr_err_vec{k}(1)]; %added by macshen
    [estusr{k},H{k}] = olspos(prvec{k},svxyzmat{k});  %sovle for the vehicle pose given pr measurement, macshen
    
    estenu{k}(i,:) = ( xyz2enu(estusr{k}(1:3),orgxyz) )';   %transform the ECEF coordinate to local one with origin orgxyz
    err{k}(i,1:3) = estenu{k}(i,1:3) - usrenu{k}(i,:);
    terr{k}(i) = estusr{k}(4);  % true clk bias is zero
    end   %end for loop of vehicle #
    
%     P0=inv(H{1}.'*H{1});
%     pf=zeros(Np,3); %first two colomns are coordinates, 3rd colomn is weight,
%     %here particles represents the offset vector that correct the vehicles' pose, therefore the particles are common to each vehicle
%     %particles are initilized every time step
%     pf(:,3)=1/Np;   %initialize weight
%     pf(:,1:2)=mvnrnd([0,0],9*P0(1:2,1:2),Np);   %try Gaussian sampling, this should be modified later


    
    
    figure(1);
    clf;
    hold on
    
%    plotcov2d(estenu{1}(i,1),estenu{1}(i,2),[1,0;0,1],'b',0,0,0,1);  %plot the GPS estimated position of vehicle as a circle
%    plotSamples([pf(:,1)+estenu{1}(i,1),pf(:,2)+estenu{1}(i,2)],'b');  %plot the hypothetical position of 1st vehicle before particle filter
    
    for k=1:N  %pf loop for each vehicle
        sigma_multipath2=1.67*noise_scaling_factor(4)^2;   %variance of mutlipath,calculated from the multipath noise in the simulation
        sigma_thermal2=1.0*noise_scaling_factor(1)^2;
        sigma_eta2=sigma_multipath2+sigma_thermal2;
        sigma_map2=0.01;  %variance of map inprecision, this might be larger for real map
        
        sigma_x=0.1*deltat^4*0.5^2;  %unmodeled position noise caused by acceleration,assuming the variance of acceleration is 1.
        sigma_x_dot=0.1*deltat^2;  %unmodeled velocity noise caused by acceleration;
        
if i==1   %for the first time step, initialize particles for ego-state estimation and multipath mitigation
    initialize_pf_ego(Np_pose,[usrenu{1}(1,1),usrenu{1}(1,2),0],[velocity,0,0],2);
    pf_update_ego(prvec{k},svxyzmat{k},sigma_thermal2,orgxyz);   %pf update given measurement
else  %second step and so on
    sigma_b=0.001*deltat^4*0.5^2;
    sigma_b_dot=0.001*deltat^2;   %specify the uncertainty level of the clock bias
    pf_prediction_ego(sigma_x,sigma_x_dot,sigma_b,sigma_b_dot,deltat);  %in this step, the multipath is not changed, it should be considered whether this will cause over confidence.
    pf_update_ego(prvec{k},svxyzmat{k},sigma_thermal2,orgxyz);
end
%after calculating weight, resample
    resample_ego;
variance(i)=pf_pose.mpvar(1,1);

%         P=inv(H{k}.'*H{k})*sigma_eta2;  %whole position covariance caused by non-common noise
%         P=P(1:2,1:2);  %extract the sub covariance matrix of x,y position
%         
%         if i==1  %initialize estimation using first measurement
%         state_mean{k}(1:2,1)=estenu{k}(i,1:2)';
%         state_mean{k}(3:4,1)=[0;0];  %set initial velocity as zero
%         state_cov{k}=eye(4,4);    %initialize state vector and its covariance matrix for the Kalman filter
%         state_cov{k}(1:2,1:2)=P(1:2,1:2);
%         else %for second and later loop, use Kalman filter
%         [state_mean{k},state_cov{k}]=Kalman_prediction(state_mean{k},state_cov{k},sigma_x,sigma_x_dot,deltat);  
%         [state_mean{k},state_cov{k}]=Kalman_update(state_mean{k},state_cov{k},estenu{k}(i,1:2)',P);
%         end
%         
%         
%         pf(:,1:2)=pf(:,1:2)+0.2*randn(Np,2);  %This line adds some random noise to the particles to increase diversity
%         %while it should be noted that this would cause a wrong posterior
%         for j=1:Np %loop for particles' weights
%         %pf(j,3)=cal_weight(pf(j,1:3),estenu{k}(i,1:2),P,sigma_map2,lane_width,vehicle_width,grid_size); %the last argument allow noise to be added to the particles
%         pf(j,3)=cal_weight(pf(j,1:3),state_mean{k}(1:2,1)',state_cov{k}(1:2,1:2),sigma_map2,lane_width,vehicle_width,grid_size); %calculate weight
%         end
%         pf=resample(pf);  %resample, %weight needn't be normalized. After resample, weight has been normalized
    end   %end (for k=1:N)
    
%     [mu,cov]=pf_to_Gaussian(pf);  %transform particle filter to Gaussian
%     deter(i)=det(cov);  %calculate the determinant of the covariance matrix
%     
%     err_temp=[];
%     for k=1:N
%     err_temp=[err_temp;err{k}(i,1:2)];   %extract the error vector for each vehicle
%     end
%     deter_error(i)=cal_deter_error(err_temp);  %calculate the determinant of the covariance of error
    
    %calculate the determinant of the variance of the true error, this is expected to be related with the estimation covariance
    plotcov2d(usrenu{1}(i,1),usrenu{1}(i,2),[1,0;0,1],'r',0,0,0,1);  %plot the true position of vehicle as a circle
    for m=1:pf_pose.N
    plotSamples(pf_pose.state(m,1:2),'b');   %plot particle for ego localization
    end
    %when not Kalman
    %plotSamples([pf(:,1)+estenu{1}(i,1),pf(:,2)+estenu{1}(i,2)],'r');  %plot the estimated position sample of 1st vehicle
    %plotcov2d(mu(1)+estenu{1}(i,1),mu(2)+estenu{1}(i,2),cov,'g',0,0,0,3);  %plot the true position of vehicle as a circle
    %when use Kalman filter
    %plotSamples([pf(:,1)+state_mean{1}(1,1),pf(:,2)+state_mean{1}(2,1)],'r');  %plot the estimated position sample of 1st vehicle
    %plotcov2d(mu(1)+state_mean{1}(1,1),mu(2)+state_mean{1}(2,1),cov,'g',0,0,0,3);  %plot the true position of vehicle as a circle
    axis equal;
    pause(1);
    hold off
 	% waitbar(i/EndLoop)   
end  %end (for i = 1:Ns)
%close(bar1);
% x=100*rand(1000,1);
% y=100*rand(1000,1);
% for k=1:1000
% d(k)=Distance_to_road([x(k),y(k)],grid_size);    %pay attention to changing the grid-size if 
% end