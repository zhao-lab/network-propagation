%  sndemo40.m     Point Positioning with RINEX2 ephemeris and observation files
%                 C/A-Code pseudoranges are used; Broadcast IONO correction and
%                 standard TROPO correction applied
%
%	M. & S. Braasch 08-2003; 07-2006
%	Copyright (c) 2003; 2006 by GPSoft
%	All Rights Reserved.
%
clear all
close all
%
global SQRTSMA
global SVID_MAT TOWSEC PHASE1 PHASE2 C1 P1 P2 D1 D2
global PHASE1LLI PHASE1SS PHASE2LLI PHASE2SS
global C1LLI C1SS P1LLI P1SS P2LLI P2SS
global MARKER_XYZ ANTDELTA OBSINT TIMESTART TIMESTOP CLOCKOFFSET
%
% load day101snav
% load mthn101sob
% 
% mthnxyz = [901906.698 -5726170.573 2651517.227];   % Marathon, Florida  CORS station

%   load 924_ding_nav
%   load 924_ding_obs
%  load multipath_ding_nav
%  load multipath_ding_obs
 load hp_1014_nav
 load hp_1014_obs


 %first eliminate the data for G04 since the ephemeris is unavailable
 N=size(C1,2);  %get the length of the data
 Z=zeros(1,N);
 C1(4,:)=Z;
 C1LLI(4,:)=Z;
 C1SS(4,:)=Z;
 D1(4,:)=Z;
 PHASE1(4,:)=Z;
 PHASE1LLI(4,:)=Z;
 PHASE1SS(4,:)=Z;
 S1(4,:)=Z;
 SVID_MAT(4,:)=Z;
 
 C1(7,:)=Z;
 C1LLI(7,:)=Z;
 C1SS(7,:)=Z;
 D1(7,:)=Z;
 PHASE1(7,:)=Z;
 PHASE1LLI(7,:)=Z;
 PHASE1SS(7,:)=Z;
 S1(7,:)=Z;
 SVID_MAT(7,:)=Z;
 
  C1(26,:)=Z;
 C1LLI(26,:)=Z;
 C1SS(26,:)=Z;
 D1(26,:)=Z;
 PHASE1(26,:)=Z;
 PHASE1LLI(26,:)=Z;
 PHASE1SS(26,:)=Z;
 S1(26,:)=Z;
 SVID_MAT(26,:)=Z;
 
   C1(16,:)=Z;
 C1LLI(16,:)=Z;
 C1SS(16,:)=Z;
 D1(16,:)=Z;
 PHASE1(16,:)=Z;
 PHASE1LLI(16,:)=Z;
 PHASE1SS(16,:)=Z;
 S1(16,:)=Z;
 SVID_MAT(16,:)=Z;
%  
%  SVID_MAT(2,:)=Z;
%  SVID_MAT(30,:)=Z;
%  SVID_MAT(21,:)=Z;
%  SVID_MAT(18,:)=Z;
 
%     C1(2,:)=Z;
%  C1LLI(2,:)=Z;
%  C1SS(2,:)=Z;
%  D1(2,:)=Z;
%  PHASE1(2,:)=Z;
%  PHASE1LLI(2,:)=Z;
%  PHASE1SS(2,:)=Z;
%  S1(2,:)=Z;
%  SVID_MAT(2,:)=Z;
%  
%     C1(30,:)=Z;
%  C1LLI(30,:)=Z;
%  C1SS(30,:)=Z;
%  D1(30,:)=Z;
%  PHASE1(30,:)=Z;
%  PHASE1LLI(30,:)=Z;
%  PHASE1SS(30,:)=Z;
%  S1(30,:)=Z;
%  SVID_MAT(30,:)=Z;
 
%mthnxyz = [901906.698 -5726170.573 2651517.227];   % Marathon, Florida  CORS station
%usrxyz = mthnxyz;
usrxyz=llh2xyz([42.29728877393/180*pi,-83.70654541167/180*pi,266.6204]);
%%RTK position
%usrxyz=llh2xyz([42.29947006458/180*pi,-83.70547232203/180*pi,270.6485]);
%usrxyz=llh2xyz([42.2971055/180*pi,-83.706585/180*pi,265.0]);

usr_clk = 0;
sv_clk = zeros(1,32);
c = 299792458;
maskangle = 5;

[r,num_epochs] = size(C1);
bar1 = waitbar(0,'Calculating Position...   ');
n = 0;

% %added by macshen to store the time history of pr
% pr_stack=[];
m=0;
for i = 1:1:num_epochs,  %Only processing every 5th sample to speed up this demo
    n = n + 1;   % Counter for the number of samples being processed
    time = TOWSEC(i);     % Time of reception given in GPS time-of-week in seconds
    clear id prvec svxyzmat
    id = find(SVID_MAT(:,i)==1);
    k = 0;
    for j = 1:length(id),   % Loop over all satellites being tracked during this epoch
        if SQRTSMA(id(j)) < 1,
            error('Ephemeris not available for all satellites being tracked')
        end
      %  if (C1(id(j),i)>SQRTSMA(id(j)))&(P1(id(j),i)>SQRTSMA(id(j)))&(P2(id(j),i)>SQRTSMA(id(j))),     % skip over obviously bad measurements
            tot_est = time -(1/c)*(C1(id(j),i)-usr_clk+sv_clk(id(j)));   % Estimate the time-of-transmission
            [svxyz,E] = svposeph(id(j),tot_est);      % Calculate the position of the satellite
            [sv_clk(id(j)),grpdel] = svclkcorr(id(j),tot_est,E);    % Calculate the satellite clock correction
            if i == 1, estusr(1:3) = usrxyz; end
            svenu = xyz2enu(svxyz,estusr(1:3));     % Convert the satellite position to east-north-up (i.e., local level) coordinates
            el = (180/pi)*atan(svenu(3)/norm(svenu(1:2)));
            if el >= maskangle,
               k = k + 1;       % counter for the satellites being utilized in this epoch
               prCA = C1(id(j),i) + sv_clk(id(j)) - grpdel*c;   % C/A-code pseudorange corrected for satellite clock and Tgd
               prvec(k) = prCA;
               if i > 1,
                   ionod = ionocorr(tot_est,svxyz,estusr(1:3));   % Calculate the broadcast IONO model correction
                   tropd = tropocorr(svxyz,estusr(1:3));    % Calculate the tropospheric correction
                   prvec(k) = prvec(k) - ionod - tropd;     % Adjust the pseudorange for the iono and tropo corrections
               end
               svxyzr = erotcorr(svxyz,prvec(k));   % Adjust satellite position coordinates for earth rotation correction
               svxyzmat(k,:) = svxyzr';
               svvis(id(j),i) = 1;
            end
     %   end

     
     
     if id(j)==21
         true_pr21(i)=norm(svxyzr'-usrxyz);
         pr21(i)=prvec(k);
     end
          if id(j)==29
         true_pr29(i)=norm(svxyzr'-usrxyz);
         pr29(i)=prvec(k);
          end
               if id(j)==20
         true_pr20(i)=norm(svxyzr'-usrxyz);
         pr20(i)=prvec(k);
               end
               if id(j)==18
         true_pr18(i)=norm(svxyzr'-usrxyz);
         pr18(i)=prvec(k);
               end
               if id(j)==15
         true_pr15(i)=norm(svxyzr'-usrxyz);
         pr15(i)=prvec(k);
               end
               if id(j)==13
         true_pr13(i)=norm(svxyzr'-usrxyz);
         pr13(i)=prvec(k);
               end
               if id(j)==5
         true_pr5(i)=norm(svxyzr'-usrxyz);
         pr5(i)=prvec(k);
               end
                    if id(j)==2
         true_pr2(i)=norm(svxyzr'-usrxyz);
         pr2(i)=prvec(k);
                    end
                         if id(j)==30
         true_pr30(i)=norm(svxyzr'-usrxyz);
         pr30(i)=prvec(k);
     end
     
    end
    % if length(prvec)==5
        m=m+1;
     svid4_all{m}=id; 
     svxyz4_all{m}=svxyzmat;
     prvec4_all{m}=prvec.';

    %end
%     %added by macshen to calculate the true range based on RTK measurement
%     true_pr(k,i)=norm(svxyzmat(k,:)-usrxyz);
%     pr_stack=[pr_stack,prvec.'];
    %
    if length(prvec) < 4,
        enuerr(n,:) = NaN;
    else
        estusr = olspos(prvec,svxyzmat);    % Ordinary least-squares position solution
        enuerr(n,:) = ( xyz2enu_exp(estusr(1:3)) )';
        usr_clk = estusr(4);
        usr_clk_vec(n) = usr_clk;
    end
    timeofsamples(n) = time;
    %
    time2_all(m)=time;
    pose2_all(:,m)=enuerr(n,1:2)';
    waitbar(i/num_epochs,bar1)
end
close(bar1)
close all
plot(enuerr(:,1)-1.98,enuerr(:,2)-5.2582,'*')
axis('square')
axis('equal')
axis([-10 10 -10 10])
grid
title('GPS C/A-Code Point Positioning Error')
ylabel('north error (m)')
xlabel('east error (m)')
text(-7,-3,'TROPO and Broadcast IONO corrections applied')

dif21=pr21-true_pr21;
dif29=pr29-true_pr29;
dif=dif21-dif29;
err_norm=sqrt(enuerr(:,1).^2+enuerr(:,2).^2);