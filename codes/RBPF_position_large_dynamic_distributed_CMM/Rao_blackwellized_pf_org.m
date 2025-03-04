% This is a test program for Rao Blackwellized particle filter
% for CMM localization. The methodology is to use particle to represent the
% distribution of the common pseudorange error while use EKF to represent
% the conditional vehicle pose distribution. Since the vehicle poses are
% conditionally independent given the common error. The dimension of the
% EKF is relatively small.
% This filter is actually centralized since the particle filter maintains
% all the information available. The map matching is treated as measurement
% to determine the weight of each particle.
clear all
load medium_sparse.mat;

file_path = './support_files';
addpath(file_path)

for tt=1:1
    %created by macshen
    tt
    clearvars -except record_err tt record_cov new_vehicle
    
    N = length(new_vehicle(:,1));
    Radius = 3000;
    for k=1:N
        for j=1:N
            d_vehicle(k,j) = norm(new_vehicle(k,1:2)-new_vehicle(j,1:2));
        end
    end
    A_con = d_vehicle < Radius;
     
    close all
    global pf;
    %N=50; % # of total vehicles
    map_angle = new_vehicle(:,3);    % randomly generate map_angle
    
    Ng = 2000;   % # of discrete grids of multipath estimation for each lane
    mp_gridsize = 0.25;   % gridsize*Ng=lane length
    lane_width = 3.5;   % width of a single lane
    vehicle_width=1.8;  % use vehicle width to increase the positioning accuracy, it is assumed that the GPS receiver locates at the center of the vehicle
    velocity=0.01;  % vehicle velocity
    Ns=100;  % # of simulation time points
    Nsv=6;  %# of visible satellites
    Np=50; % # of particles
    block_prob=0;
    
    %A_con=(eye(N)+rand(N))>0.75;
    
    %redundant
    mpmat=cell(N,1);   %multipath error time history for N vehicles
    for i=1:N
        %generate multipath error, use different random seeds
        mpmat{i}=mpgen(24,1200,2,floor(10000*rand(1))+54321+i);
    end
    %mpmat{1}(1:Ns,1)=10000*sin(1/100*(1:Ns))';
    % Origin of the local coordinate system, first latitude, second longitudinate:
    orgllh = [40*pi/180 80*pi/180 0];
    % Convert the origin to the xyz ECEF cooridinate system
    orgxyz = llh2xyz(orgllh);
    loadgps
    startt=500; t = startt; deltat=0.1;
    
    % usrenu{1}(1:Ns,1)=(6:velocity*deltat:6+velocity*deltat*(Ns-1))';
    % usrenu{1}(1:Ns,1)=0.5*((1:Ns)*deltat).^2;
    %
    % usrenu{1}(1:Ns,2)=50-lane_width/2;
    % usrenu{1}(1:Ns,3)=0;
    for k=1:N
        usrenu{k}(1:Ns,1:3)=0;
    end
    % for k=1:1
    % usrenu{k}(1:Ns,1)=(16:velocity*deltat:16+velocity*deltat*(Ns-1))';
    % usrenu{k}(1:Ns,2)=250-lane_width/2;
    % usrenu{k}(1:Ns,3)=0;
    % end
    % for k=2:2
    % usrenu{k}(1:Ns,1)=(473:-velocity*deltat:473-velocity*deltat*(Ns-1))';
    % usrenu{k}(1:Ns,2)=250+lane_width/2;
    % usrenu{k}(1:Ns,3)=0;
    % end
    % for k=3:3
    % usrenu{k}(1:Ns,2)=(16:velocity*deltat:16+velocity*deltat*(Ns-1))';
    % usrenu{k}(1:Ns,1)=250+lane_width/2;
    % usrenu{k}(1:Ns,3)=0;
    % end
    % for k=4:4
    % usrenu{k}(1:Ns,2)=(473:-velocity*deltat:473-velocity*deltat*(Ns-1))';
    % usrenu{k}(1:Ns,1)=250-lane_width/2;
    % usrenu{k}(1:Ns,3)=0;             %generate the vehicles' path in the local coordinate system
    % end
    %
    % for k=5:N
    %     index=mod(k,4);
    %     if index==0
    %         index=4;
    %     end
    %     usrenu{k}=usrenu{index};
    % end
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
    
    
    % load('mp_x.mat');
    % load('mp_y.mat');
    load('mp_30.mat');
    global mp_xg;
    global mp_yg;
    global LOS_xg;
    global LOS_yg;   %global variable shared with add_mp_wn
    
    for k=1:Nsv
        mp_x{k}=mp_x{k}+0.1*abs(randn(size(mp_x{k})));
        mp_y{k}=mp_y{k}+0.1*abs(randn(size(mp_y{k})));   %perturb to simulate diffraction
        wavelength=3*10^8/1575.42/10^6;   % L1 signal wavelength
        phase_x=2*pi*mp_x{k}/wavelength;
        phase_y=2*pi*mp_y{k}/wavelength;
        mp_x{k}=0.5*mp_x{k}.*cos(phase_x);
        mp_y{k}=0.5*mp_y{k}.*cos(phase_y);
    end
    mp_xg=mp_x;
    mp_yg=mp_y;
    LOS_xg=LOS_x;
    LOS_yg=LOS_y;
    for k=1:6
        mp_xg{k}=zeros(1920,40);
        mp_yg{k}=zeros(1920,40);
        LOS_xg{k}=ones(1920,40);
        LOS_yg{k}=ones(1920,40);
    end
    
    for k=1:N
        com_list{k}=k;
    end
    
    for i = 1:Ns  %draw the first step
        
        %     if mod(i,1)==0   %this is not the fusion frequency
        %     for k=1:N
        %         while 1
        %             new=ceil(rand(1)*N);   %generate a new com that is different from k
        %             if new~=k
        %                 break;
        %             end
        %         end
        % %         if k==1    %prevent the 1st vehicle connecting to the 4th one
        % %             if rand(1)<0.5
        % %                 new=2;
        % %             else
        % %                 new=3;
        % %             end
        % %         end
        %     new_list{k}=new;
        %     end
        %     end
        
        %     for k=1:N
        %         new_added{k}=0;%indicate if there is a new vehicle added
        %         if ~ismember(new_list{k},com_list{k})   %so if the new_list does not change, the com_list wont be updated
        %             new_added{k}=1;
        %             if length(com_list{k})>=2
        %                 com_list{k}=[com_list{k}(1:end-1),new_list{k}];
        %             else
        %                 com_list{k}=[com_list{k},new_list{k}];
        %             end
        %         end
        %     end
        
        %     com_list{1}=[1,2];
        %     com_list{2}=[2,3];
        %     com_list{3}=[3,4];
        %     com_list{4}=[4,1];
        for k=1:N
            com_list{k}=find(A_con(k,:)==1);
            %     if length(com_list{k})>=6
            %         com_list{k}=com_list{k}(1:6);
            %     end
            if ~ismember(k,com_list{k})
                com_list{k}=[com_list{k},k];
            end
        end
        
        A_sparse=zeros(N,N);
        for k=1:N
            for j=1:length(com_list{k})
                A_sparse(k,com_list{k}(j))=1;
            end
        end      
     
        t = t + deltat;
        
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
            
            estenu{k}(i,:) = ( xyz2enu(estusr{k}(1:3),orgxyz) )';   %transform the ECEF coordinate to local one with origin orgxyz
            common_error_position(:,i)=estenu{k}(i,1:2)'-usrenu{k}(1,1:2)';
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
        
        
        % figure(1);
        clf;
        hold on
        
        
        %     for k=1:N  %pf loop for each vehicle
        sigma_multipath2=1.67*noise_scaling_factor(4)^2;   %variance of mutlipath,calculated from the multipath noise in the simulation
        %         sigma_eta2=sigma_multipath2+sigma_thermal2;
        %         sigma_map2=0.01;  %variance of map inprecision, this might be larger for real map
        
        
        if i==1   %for the first time step, initialize particles for ego-state estimation and multipath mitigation
            initialize_pf_CMM(Np,Nsv,N,common_error(i,:).',velocity);
            for j=1:N
                for v=1:length(com_list{j})
                    k=com_list{j}(v);   %for the j-th pf, update the k-th vehicle
                    if rand>block_prob
                        update_pf_CMM(estenu{k}(i,1:2)',svxyzmat{k},sigma_thermal2,orgxyz,Np,k,hor_DLP{k},H{k},j);   %pf update given measurement
                        resample_CMM(j);  %think about it, resample for every vehicle update or for the whole update?
                        weight_CMM(k,j,map_angle(k));  %calculate weight according to the map constraints
                        resample_CMM(j);
                    end
                end
            end
        else  %second step and so on
            predict_pf_CMM(sigma_a2,sigma_b2,sigma_d2,deltat,N,Np,Nsv,hor_DLP{k});
            
            
            if mod(i,1)==0
                
                for k=1:N
                    com_err(1:2,k)=pf_get_common(k,H{k});   %project the pr into common error
                end
                mean_com=mean(com_err')';
                for k=1:N
                    com_err(1:2,k)=com_err(1:2,k)-mean_com;
                end
                
                for k1=1:N
                    for k2=1:N
                        M_opt(k1,k2)=com_err(:,k1).'*com_err(:,k2);  %matrix fed to the optimization program
                    end
                end
                
                % For optimized distributed RBPF (probably with poor results; comment this if implement
                % randomly distributed RBPF
                for k=1:N
                    N_alloc{k}=opt_com_weight(M_opt,k,Np,com_list{k});   %calculate the number of particles for allocation
                end
                
                 % For random distributed RBPF (probably with poor results; comment this if implement
                 % optimated ditributed RBPF
%                 for k=1:N
%                     N_alloc{k}=round(A_con(k,:)'.*rand(N,1)*200/sum(A_con(k,:)'));
%                     N_alloc{k}(k)=Np;
%                 end
                
                N_trans_particle(i)=0;
                for k=1:N
                    N_trans_particle(i)=N_trans_particle(i)+sum(N_alloc{k});
                end
                
                for k=1:N    %merge every a fixed number of steps
                    
                    %         if length(com_list{k})>=2
                    %             for v=2:length(com_list{k})
                    %                 index=com_list{k}(v);
                    %                 pf{k}=[pf{k},pf{index}];
                    %             end
                    %         end
                    
                    for j=1:N
                        if j~=k
                            if N_alloc{k}(j)>0
                                numb=N_alloc{k}(j);    %numb might be larger than the total number of particle
                                while numb>Np
                                    pf{k}=[pf{k},pf{j}];
                                    numb=numb-Np;
                                end
                                pf{k}=[pf{k},pf{j}(1:numb)];
                            end
                        end
                    end
                    
                end
                
                for k=1:N
                    for j=1:N
                        weight_CMM(j,k,map_angle(j));
                    end
                    resample_CMM_with_desample(k,Np);  %here need consideration!!!!!!!!!!!
                end
                
            end
            
            
            for j=1:N
                for v=1:length(com_list{j})
                    k=com_list{j}(v);   %for the j-th pf, update the k-th vehicle
                    if rand>block_prob
                        update_pf_CMM(estenu{k}(i,1:2)',svxyzmat{k},sigma_thermal2,orgxyz,Np,k,hor_DLP{k},H{k},j);   %pf update given measurement
                        resample_CMM(j);
                        weight_CMM(k,j,map_angle(k));  %calculate weight according to the map constraints
                        resample_CMM(j);
                    end
                end
            end
            
            
            % for k=1:N    %merge upon new connection
            %     if new_added{k}==1
            %         index=new_list{k};
            %         pf{k}=[pf{k},pf{index}];
            %         for j=1:N
            %         weight_CMM(j,k);
            %         end
            %         resample_CMM_with_desample(k,Np);  %here need consideration!!!!!!!!!!!
            %      end
            % end
            
            
            
        end
        
        for k=1:2
            [mu{k},cov{k}]=pf_to_gaussian_RBPF(k);
        end
        
        %         pf=resample(pf);  %resample, %weight needn't be normalized. After resample, weight has been normalized
        
        %calculate the determinant of the variance of the true error, this is expected to be related with the estimation covariance
        % REMEMBER: commented only for faster simulation
        plotcov2d(usrenu{1}(i,1),usrenu{1}(i,2),[1,0;0,1],'r',0,0,0,1);  %plot the true position of vehicle as a circle
        for k=1:Np
            plotcov2d(pf{1}(k).mu{1}(1),pf{1}(k).mu{1}(3),pf{1}(k).cov{1}([1,3],[1,3]),'b',0,0,0,3);  %plot the estimated position of vehicle as a circle
        end
        plotcov2d(mu{1}(1),mu{1}(2),cov{1},'g',0,0,0,3);
        
        ave_err_CMM(i,1:2)=zeros(1,2);
        var_err_CMM(i)=0;
        ave_sq_err(i)=0;
        for k=1:N
            err_CMM{k}(i,:)=mu{k}'-usrenu{k}(i,1:2);
            ave_err_CMM(i,1:2)=ave_err_CMM(i,1:2)+err_CMM{k}(i,:)/N;
            bias_x{k}(i)=err_CMM{k}(i,1);
            bias_y{k}(i)=err_CMM{k}(i,2);
            err_norm{k}(i)=norm(err_CMM{k}(i,1:2));
            deter{k}(i)=det(cov{k});
        end
        for k=1:N
            var_err_CMM(i)=var_err_CMM(i)+norm((err_CMM{k}(i,:)-ave_err_CMM(i,1:2))')^2/N;
            ave_sq_err(i)=ave_sq_err(i)+norm(err_CMM{k}(i,:)')^2/N;
            norm_ave_err(i)=norm(ave_err_CMM(i,1:2)');
        end
        rt_var_err(i)=sqrt(var_err_CMM(i));
        rt_ave_sq_err(i)=sqrt(ave_sq_err(i));
        
        %when not Kalman
        %plotSamples([pf(:,1)+estenu{1}(i,1),pf(:,2)+estenu{1}(i,2)],'r');  %plot the estimated position sample of 1st vehicle
        %plotcov2d(mu(1)+estenu{1}(i,1),mu(2)+estenu{1}(i,2),cov,'g',0,0,0,3);  %plot the true position of vehicle as a circle
        %when use Kalman filter
        %plotSamples([pf(:,1)+state_mean{1}(1,1),pf(:,2)+state_mean{1}(2,1)],'r');  %plot the estimated position sample of 1st vehicle
        %plotcov2d(mu(1)+state_mean{1}(1,1),mu(2)+state_mean{1}(2,1),cov,'g',0,0,0,3);  %plot the true position of vehicle as a circle
        axis equal;
        pause(0.0001);
        % waitbar(i/EndLoop)
    end  %end (for i = 1:Ns)
    %close(bar1);
    % x=100*rand(1000,1);
    % y=100*rand(1000,1);
    % for k=1:1000
    % d(k)=Distance_to_road([x(k),y(k)],grid_size);    %pay attention to changing the grid-size if
    % end
    %
    for i = 1:N
        record_err(tt,i)=mean(err_norm{i});
        record_cov(tt,i)=mean(deter{i});
    end
end

%% PLOT ALL ERRORS AND CORVARIANCE
figure;
plot(rt_var_err,'linewidth',2)
title('Different types of CMM error')
hold on;
plot(rt_ave_sq_err,'linewidth',2)
hold on;
plot(var_err_CMM,'linewidth',2)
hold on;
plot(ave_sq_err,'linewidth',2)
hold on;
plot(norm_ave_err,'linewidth',2)
legend('rt-var-err','rt-ave-sq-err','var-err-CMM','ave-sq-err','norm-ave-err')
figure;
plot(ave_err_CMM,'linewidth',2)
legend('ave-err-CMM_x','ave-err-CMM_y')
title('Average Error for CMM_x and CMM_y')
% figure;
% plot(common_error,'linewidth',2)
% legend('common-error_1','common-error_2','common-error_3','common-error_4', ...
%     'common-error_5','common-error_6')
figure;
plot(common_error_position','linewidth',2)
legend('common-error-position_x','common-error-position_y')
title('Common Error for Longitude and Lateral Position')
figure;
plot(com_err(1,:),'-o','linewidth',2);
hold on;
plot(com_err(2,:),'-o','linewidth',2);
legend('com-err_x','com-err_y')
title('Com-Err for cars in network')
%%
for i = 1:N
    for j = 1:Ns
         err_CMM_norm{i}(j) = norm(err_CMM{i}(j,1),err_CMM{i}(j,2));
    end
end
% figure;
 for i = 1:10:N
    figure;
    plot(err_CMM_norm{i},'linewidth',2);
%     hold on;
 end
