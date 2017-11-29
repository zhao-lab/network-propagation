%this is a test program for Rao Blackwellized particle filter
%for CMM localization. The methodology is to use particle to represent the
%distribution of the common pseudorange error while use EKF to represent the
%conditional vehicle pose distribution. Since the vehicle poses are conditionally
%independent given the common error. The dimension of the EKF is relatively small.
%This filter is actually centralized since the particle filter maintains
%all the information available. The map matching is treated as measurement
%to determine the weight of each particle.
clear all
load vehicle_50.mat;
% load most_sparse.mat
% vehicle=new_vehicle;

file_path = './support_files';
addpath(file_path)

% if medium or most sparse, add the followinng
% vehicle = new_vehicle;

for tt=1:1
    %created by macshen
    tt
    clearvars -except record_err tt vehicle record_err_CMM
    close all
    global pf;
    map_angle = vehicle(:,3);    %randomly generate map_angle
    N = length(vehicle(:,1)); % # of total vehicles
    % load('distance_and_gradient.mat');  %contains grid_size, distance and gradient of distance from road centers
    % global dis;
    % global grad;   %global quantities to share with sub-function for interpolating distances
    % dis=distance;
    % grad=grad_dis;
    
    Ng=2000;   % # of discrete grids of multipath estimation for each lane
    mp_gridsize=0.25;   % gridsize*Ng = lane length
    lane_width=3.5;   %width of a single lane
    vehicle_width=1.8;  %use vehicle width to increase the positioning accuracy, it is assumed that the GPS receiver locates at the center of the vehicle
    %     velocity=100*randn;  %vehicle velocity
    Ns=100;  % #of simulation time points
    Nsv=6;  %#of visible satellites
    Np=50; % #of particles
    %     block_prob=0;
    % randomly block a small amount of neighboring vehicles;
    block_prob=0;
    
    %define a dynamic network system with changing distance wrt time and
    %velocity
    velocity_norm = 10*randn(1,N);
    
    velocity_norm(1) = 0;
    velocity_norm(2) = 5;
    velocity_norm(3) = 20;
    
    velocity = zeros(2,N);
    for i = 1:N
        velocity(1,i) = velocity_norm(i)*cos(vehicle(i,3)); %0;
        velocity(2,i) = velocity_norm(i)*sin(vehicle(i,3)); %0;
    end
    
    distance = zeros(N,N);
    %     distance_dyn = zeros(N,N,Ns);
    for k = 1:Ns
        for i = 1:N
            for j = 1:N
                distance(i,j) = sqrt((vehicle(i,1)-vehicle(j,1))^2 ...
                    + (vehicle(i,2) - vehicle(j,2))^2);
                %                 distance_dyn(i,j,k) = distance(i,j) + ...
                %                     (velocity(i)-velocity(j))*Ns;
                distance_dyn{k}(i,j) = sqrt((vehicle(i,1)-vehicle(j,1)+(velocity(1,i)-velocity(1,j))*1*k)^2 ...
                    + (vehicle(i,2) - vehicle(j,2)+(velocity(2,i)-velocity(2,j))*1*k)^2);
                position{k} = vehicle(:,1:2);
                position_dyn{k}(i,1) = position{k}(i,1) + velocity(1,i)*k;
                position_dyn{k}(i,2) = position{k}(i,2) + velocity(2,i)*k;
            end
        end
    end
    %%
    %redundent
    mpmat=cell(N,1);   %multipath error time history for N vehicles
    for i=1:N
        mpmat{i}=mpgen(24,1200,2,floor(10000*rand(1))+54321+i); %generate multipath error, use different random seeds
    end
    %mpmat{1}(1:Ns,1)=10000*sin(1/100*(1:Ns))';
    orgllh = [40*pi/180 80*pi/180 0];  %origin of the local coordinate system, first latitude, second longitudinate
    orgxyz = llh2xyz(orgllh);   %convert the origin to the xyz ECEF cooridinate system
    loadgps
    startt=500; t = startt; deltat=0.1;
    
    %define usr in a complicated way
    %     usrenu{1}(1:Ns,1)=(6:velocity*deltat:6+velocity*deltat*(Ns-1))';
    %     usrenu{1}(1:Ns,1)=0.5*((1:Ns)*deltat).^2;
    %
    %     usrenu{1}(1:Ns,2)=50-lane_width/2;
    %     usrenu{1}(1:Ns,3)=0;
    %     for k=1:1
    %     usrenu{k}(1:Ns,1)=(16:velocity*deltat:16+velocity*deltat*(Ns-1))';
    %     usrenu{k}(1:Ns,2)=250-lane_width/2;
    %     usrenu{k}(1:Ns,3)=0;
    %     end
    %     for k=2:2
    %     usrenu{k}(1:Ns,1)=(473:-velocity*deltat:473-velocity*deltat*(Ns-1))';
    %     usrenu{k}(1:Ns,2)=250+lane_width/2;
    %     usrenu{k}(1:Ns,3)=0;
    %     end
    %     for k=3:3
    %     usrenu{k}(1:Ns,2)=(16:velocity*deltat:16+velocity*deltat*(Ns-1))';
    %     usrenu{k}(1:Ns,1)=250+lane_width/2;
    %     usrenu{k}(1:Ns,3)=0;
    %     end
    %     for k=4:4
    %     usrenu{k}(1:Ns,2)=(473:-velocity*deltat:473-velocity*deltat*(Ns-1))';
    %     usrenu{k}(1:Ns,1)=250-lane_width/2;
    %     usrenu{k}(1:Ns,3)=0;             %generate the vehicles' path in the local coordinate system
    %     end
    %
    %     for k=5:N
    %         index=mod(k,4);
    %         if index==0
    %             index=4;
    %         end
    %         usrenu{k}=usrenu{index};
    %     end
    %end usr difinition/complex (not workable)
    
    %define usr in a simple&workable way
    for k=1:N
        usrenu{k}(1:Ns,1:3)=0;
    end
    %end usr definition/simple
    
    %EndLoop = max(size(usrenu));
    %bar1 = waitbar(0,'Generating User Solutions...  ');
    %id=[]; %added by macshen
    for k=1:N
        pr_err{k}=[]; %added by macshen
        mp_err{k}=[];
    end
    state=[];
    
    sigma_a2=1;  %unmodeled position noise caused by acceleration,assuming the variance of acceleration is 1.
    sigma_b2=1;
    sigma_d2=1;   %specify the uncertainty level of the clock bias
    noise_scaling_factor=[1 1 0 0 1];
    for k=1:N
        white_noise{k}=noise_scaling_factor(1)^2*randn(Ns,Nsv);
    end
    
    %     map_mp=mp_generation(Nsv,[0,0,0,1,0,1],Ng,mp_gridsize);%generate mp according to geometry
    %     map_mp=abs(map_mp);
    %load('mp.mat');
    map_mp=[];
    %         sigma_m2=1;   %mp variance
    %         sigma_md2=1;   %mp drift variance
    
    %     % Comment this part if multipath error is neglected
    %     load('mp_x.mat');
    %     load('mp_y.mat');
    %     load('mp_30.mat');
    %     global mp_xg;
    %     global mp_yg;
    %     global LOS_xg;
    %     global LOS_yg;   %global variable shared with add_mp_wn
    %
    %
    %     for k=1:Nsv
    %     mp_x{k}=mp_x{k}+0.1*abs(randn(size(mp_x{k})));
    %     mp_y{k}=mp_y{k}+0.1*abs(randn(size(mp_y{k})));   %perturb to simulate diffraction
    %     wavelength=3*10^8/1575.42/10^6;   % L1 signal wavelength
    %     phase_x=2*pi*mp_x{k}/wavelength;
    %     phase_y=2*pi*mp_y{k}/wavelength;
    %     mp_x{k}=0.5*mp_x{k}.*cos(phase_x);
    %     mp_y{k}=0.5*mp_y{k}.*cos(phase_y);
    %     end
    %     mp_xg=mp_x;
    %     mp_yg=mp_y;
    %     LOS_xg=LOS_x;
    %     LOS_yg=LOS_y;
    %     for k=1:6
    %         mp_xg{k}=zeros(1920,40);
    %         mp_yg{k}=zeros(1920,40);
    %         LOS_xg{k}=ones(1920,40);
    %         LOS_yg{k}=ones(1920,40);
    %     end
    
    
    %% attempt to add "range of confidence"
    
    confi_rang = 3000;
    for j = 1:Ns
        ni = 1;
        for k = 1:3
        for i = 1:N
            if distance_dyn{1,j}(i,k) < confi_rang
                dis_in_ind{j,1}(ni,k) = i;
                dis_in{j,1}(ni,k) = distance_dyn{1,j}(i,k);
                ni = ni + 1;
            end
        end
        end
    end
    % drop the ";" if you want to see how many vehicles are within the distance
    % range of higher "confidence"
    dis_in_ind;
    
    %%
    
    num_block_drop = zeros(Ns,N,3);
    block_thrend = rand(Ns,N,3);
    for i = 1:Ns  %draw the first step
        t = t + deltat;
        
        for j = 1:3
        for k=1:N   %GPS measurement loop for all the vehicles
            usrxyz{k}=enu2xyz(usrenu{k}(i,:),orgxyz);
            [svxyzmat{k},svid{k},elevation_angle{k}] = gensv(usrxyz{k},t,5,[],Nsv);   % here 5 is the mask angle to determine which satellites are visible to the vehicles, macshen
            % id=[id;svid(1)]; The last argument is to specify the wanted number of satellite and give the first group of  available satellites %added by macshen
            [prvec{k},adrvec{k},pr_err_vec{k},common_error(i,:)] = genrng_common_noise(k,usrxyz{k},svxyzmat{k},svid{k},t,noise_scaling_factor,[],mpmat{k});  %pr_err_vec added by macshen, receiver identification k seems to be redundant here. the magnitude of noise flag is thermal,Tropospheric,SA,Multipath,Ionospheric
            [prvec{k}]=add_mp_wn(usrenu{k}(i,:),prvec{k},pr_err_vec{k},map_mp,mp_gridsize,Ng,white_noise{k}(i,:));
            [estusr{k},H{k}] = olspos_clone(prvec{k},svxyzmat{k},k,orgxyz);  %sovle for the vehicle pose given pr measurement, macshen
            %    common_error_position(:,i)=inv(H{1}'*H{1})*H{1}'*common_error(i,:)';
            hor_DLP{k}=inv(H{k}'*H{k});
            hor_DLP{k}=hor_DLP{k}(1:2,1:2);
            
            estenu{k}(i,:) = (xyz2enu(estusr{k}(1:3),orgxyz))';   %transform the ECEF coordinate to local one with origin orgxyz
            
            common_error_position_1(:,i)=estenu{k}(i,1:2)'-usrenu{k}(1,1:2)';
            common_error_position_2(:,i)=estenu{k}(i,1:2)'-usrenu{k}(2,1:2)';
            
            
            %    err{k}(i,1:3) = estenu{k}(i,1:3) - usrenu{k}(i,:);
            %    terr{k}(i) = estusr{k}(4);  % true clk bias is zero
        end   %end for loop of vehicle #
        
        %     for k=1:N
        %      pr_err{k}=[pr_err{k};pr_err_vec{k}]; %added by macshen
        %      mp_err{k}=[mp_err{k};mp_err_vec{k}];
        %     end
        %     P0=inv(H{1}.'*H{1});
        %     pf=zeros(Np,3); %first two colomns are coordinates, 3rd colomn is weight,
        %     %here particles represents the offset vector that correct the vehicles' pose, therefore the particles are common to each vehicle
        %     %particles are initilized every time step
        %     pf(:,3)=1/Np;   %initialize weight
        %     pf(:,1:2)=mvnrnd([0,0],9*P0(1:2,1:2),Np);   %try Gaussian sampling, this should be modified later
        sigma_thermal2=1.0*noise_scaling_factor(1)^2;
        
        
        figure(1);
        clf;
        hold on
        
        
        %     for k=1:N  %pf loop for each vehicle
        sigma_multipath2=1.67*noise_scaling_factor(4)^2;   %variance of mutlipath,calculated from the multipath noise in the simulation
        %         sigma_eta2=sigma_multipath2+sigma_thermal2;
        %         sigma_map2=0.01;  %variance of map inprecision, this might be larger for real map
        
        
        if i==1   %for the first time step, initialize particles for ego-state estimation and multipath mitigation
            initialize_pf_CMM_new(j,Np,Nsv,N,common_error(i,:).',distance_dyn{i}(j,k));
            for k=1:N
                update_pf_CMM_new(j,estenu{k}(i,1:2)',svxyzmat{k},sigma_thermal2,orgxyz,Np,k,hor_DLP{k},H{k});   %pf update given measurement
                resample_CMM_new(j,distance_dyn{i}(j,k),confi_rang,size(dis_in{i},1),0);  %think about it, resample for every vehicle update or for the whole update?
                weight_CMM_new(j,k,map_angle(k),distance_dyn{i}(1,k));  %calculate weight according to the map constraints
                resample_CMM_new(j,distance_dyn{i}(j,k),confi_rang, size(dis_in{i},1),0);
            end
        else  %second step and so on
            predict_pf_CMM_new(j,sigma_a2,sigma_b2,sigma_d2,deltat,N,Np,Nsv,hor_DLP{k});
            for k=1:N
                if  block_thrend(i,k,j) < block_prob && distance_dyn{i}(j,k)>2000
                    %num_block = num_block+1;
                    num_block_drop(i,k) = num_block_drop(i,k) + 1;
                end
                if block_thrend(i,k,j) >= block_prob || distance_dyn{i}(j,k)<=2000
                    sum_blcok_drop(i,j) = sum(num_block_drop(i,:,j));
                    update_pf_CMM_new(j,estenu{k}(i,1:2)',svxyzmat{k},sigma_thermal2,orgxyz,Np,k,hor_DLP{k},H{k});   %pf update given measurement
                    resample_CMM_new(j,distance_dyn{i}(j,k),confi_rang,size(dis_in{i},1),sum(num_block_drop(i,:,j)));
                    %   weight_CMM(k,map_angle(k));  %calculate weight according to the map constraints
                    weight_CMM_new(j,k,map_angle(k),distance_dyn{i}(1,k));  %calculate weight according to the map constraints
                    resample_CMM_new(j,distance_dyn{i}(j,k),confi_rang,size(dis_in{i},1),sum(num_block_drop(i,:,j)));
                    % store realtime pf.weight related every timestep (i in [1,Ns])
                    % and every vehicles (k in [1,N])
                    pf_weight(i,k) = pf{j}.weight;
                end
                
            end
        end
        end
        
        
        
        %[mu,cov]=pf_to_gaussian_RBPF;
        for k=1:3
            [mu{k},cov{k}]=pf_to_gaussian_RBPF_2(k);
        end
        
        %         pf=resample(pf);  %resample, %weight needn't be normalized. After resample, weight has been normalized
        
        %%%%%%%%%%%%%%
        %     end
        %     %%
        %
        %      for i = 1:Ns
        %%%%%%%%%%%%%%
        %calculate the determinant of the variance of the true error, this is expected to be related with the estimation covariance
        plotcov2d(usrenu{1}(i,1),usrenu{1}(i,2),[1,0;0,1],'r',0,0,0,1);  %plot the true position of vehicle as a circle
        %         %% only a test:
        %         % Originally, pf(k).mu not necessarity changed
        %         for k = 1:Np
        %             pf(k).mu{1}([1,3]) = mu;
        %         end
        %         %%
        for j = 1:3
        for k=1:Np
            plotcov2d(pf{1}(k).mu{1}(1),pf{1}(k).mu{1}(3),pf{1}(k).cov{1}([1,3],[1,3]),'b',0,0,0,3);  %plot the estimated position of vehicle as a circle
            plotcov2d(pf{2}(k).mu{1}(1),pf{2}(k).mu{1}(3),pf{2}(k).cov{1}([1,3],[1,3]),'b',0,0,0,3);
            plotcov2d(pf{3}(k).mu{1}(1),pf{3}(k).mu{1}(3),pf{3}(k).cov{1}([1,3],[1,3]),'b',0,0,0,3);
        end
        % plotcov2d(mu{j}(1),mu{j}(2),cov{j},'g',0,0,0,3);
        plotcov2d(mu{1}(1),mu{1}(2),cov{1},'g',0,0,0,3);
        plotcov2d(mu{2}(1),mu{2}(2),cov{2},'k',0,0,0,3);
        plotcov2d(mu{3}(1),mu{3}(2),cov{3},'y',0,0,0,3);
        xlim([-4,4]);
        ylim([-4,4]);
        % axis equal
        err_CMM{j}(i)=norm(mu{j}'-usrenu{1}(i,1:2));
        deter{j}(i)=det(cov{j});
        %when not Kalman
        %plotSamples([pf(:,1)+estenu{1}(i,1),pf(:,2)+estenu{1}(i,2)],'r');  %plot the estimated position sample of 1st vehicle
        %plotcov2d(mu(1)+estenu{1}(i,1),mu(2)+estenu{1}(i,2),cov,'g',0,0,0,3);  %plot the true position of vehicle as a circle
        %when use Kalman filter
        %plotSamples([pf(:,1)+state_mean{1}(1,1),pf(:,2)+state_mean{1}(2,1)],'r');  %plot the estimated position sample of 1st vehicle
        %plotcov2d(mu(1)+state_mean{1}(1,1),mu(2)+state_mean{1}(2,1),cov,'g',0,0,0,3);  %plot the true position of vehicle as a circle
        % axis equal;
        pause(0.0001);
        % waitbar(i/EndLoop)
    end  %end (for i = 1:Ns)
    end
    %close(bar1);
    % x=100*rand(1000,1);
    % y=100*rand(1000,1);
    % for k=1:1000
    % d(k)=Distance_to_road([x(k),y(k)],grid_size);    %pay attention to changing the grid-size if
    % end
    for j = 1:3
    record_err{j}(tt)=mean(err_CMM{j})
    record_err_CMM{j}{tt} = err_CMM{j};
    end
end

% %% Different types of error evaluation
% % for i = 1:Ns
% %       ave_err_CMM(i,1)=zeros(1,1);
% %         var_err_CMM(i)=0;
% %         ave_sq_err(i)=0;
% %         for k=1:N
% %             err_CMM(i,:)=mu'-usrenu(i,1:2);
%             ave_err_CMM=err_CMM;
% %             bias_x(i)=err_CMM(i,1);
% %             bias_y(i)=err_CMM(i,2);
%             err_norm(i)=norm(err_CMM(i,:));
%             deter(i)=det(cov);
% %         end
%         for k=1:N
%             var_err_CMM(i)=var_err_CMM(i)+norm((err_CMM{k}(i,:)-ave_err_CMM(i,1:2))')^2/N;
%             ave_sq_err(i)=ave_sq_err(i)+norm(err_CMM{k}(i,:)')^2/N;
%             norm_ave_err(i)=norm(ave_err_CMM(i,1:2)');
%         end
%         rt_var_err(i)=sqrt(var_err_CMM(i));
%         rt_ave_sq_err(i)=sqrt(ave_sq_err(i));
% end
%% PLOT ALL ERRORS AND CORVARIANCE
% figure;
% % hold on;
% plot(common_error_position','linewidth',2)
% legend('common-error-position_x','common-error-position_y')
% % legend('common-error-position_x least sparse','common-error-position_y least sparse', ...
% %     'common-error-position_x medium sparse', 'common-error-position_y medium sparse', ...
% %     'common-error-position_x most sparse', 'common-error-position_y most sparse')
% title('Common Error for Longitude and Lateral Position')

%%
% figure;
hold on;
plot(record_err_CMM{3}{1}','.-','linewidth',2)
legend('err-CMM')
hold on;
% legend('err-CMM w/o mp')
% legend('err-CMM original', 'add dynamic and adjust update.weight w/ distance', 'add dynamic and adjust resample weight w/ distance')
% legend('error, Stationary model, equal resample weights', 'error, Dynamic model, equal resample weights(1)','error, Dynamic model, equal resample weights(2)','error, Dynamic model, filtered resample weights')
legend('w/ v_1 = 0', 'w/ v_2 = 5;','w/ v_3 = 20')
title(strcat(num2str(N), ' vehicle network, Error for CMM (v_{mean} = 23.5 m/s'))

record_err
%% Provide visualization for dynamic postions
% close all
% for each timestep, record position:
% for i = 1:100
%     figure;
%     for m = 1:N
%     % figure;
%     plot(position_dyn{i}(m,1),position_dyn{i}(m,2),'o');
%     xlabel('X/ m')
%     ylabel('Y/ m')
%     hold on;
%      % saveas(gcf,strcat('imge_', i))
%     end
% end
%%
% close all;
figure
list_1 = distance_dyn{1,1}(:,1);
hist(list_1)
title('distance distribution')
legend('timestep = 1')
%%
% hold on
figure;
list_Ns = distance_dyn{1,Ns}(:,1);
hist(list_Ns)
title('distance distribution')
legend('timestep = Ns')
%%
figure;
for i = 1:Ns
plot(i,size(dis_in_ind{i},1),'o-')
hold on;
end
dis_in_ind
%%
figure;
for i = 1:Ns
	num_dis_in(i,1) = size(dis_in_ind{i},1);
end
plot(num_dis_in,'o-')
title('number of vehicles within confident range or being droped')
hold on;
plot(sum_block_drop,'x-')
legend('number within range at each timestep','number dropped at each timestep')

dis_in_ind
sum_block_drop
%%
% close all
figure
hold on;
plot(pf_weight(Ns,:),'o-','linewidth',2)
legend('last week','modified')
title('pf-weight for particles from 24 vehicles at timestep 100')
%(1)','pf-weight from 50 vehicles (2)')
%%
%  % load the images
%  images    = cell(3,1);
%  images{1} = imread('img_16.png');
%  images{2} = imread('img_17.png');
%  images{3} = imread('img_18.png');
%  images{4} = imread('img_19.png');
%  images{5} = imread('img_20.png');
%   % create the video writer with 1 fps
%  writerObj = VideoWriter('myVideo.avi');
%  writerObj.FrameRate = 1;
%  % set the seconds per image
%  secsPerImage = [5 10 15];
%  % open the video writer
%  open(writerObj);
%  % write the frames to the video
%  for u=1:length(images)
%      % convert the image to a frame
%      frame = im2frame(images{u});
%      for v=1:secsPerImage(u)
%          writeVideo(writerObj, frame);
%      end
%  end
%  % close the writer object
%  close(writerObj);
%%

%             num_block = 0;
%             num_block_drop = zeros(N,0);
%             for k=1:N
%                 % add counter for block_prob dropped particles
%                 num_block_drop(k) = 0;
%                 if rand <= block_prob
%                     num_block = num_block+1;
%                     num_block_drop(k) = num_block_drop(k) + num_block;
%                 end
%             end

