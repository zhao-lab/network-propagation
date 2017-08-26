%  sndemo05.m       Static User   OLS positioning
clear all
close all
%    
Ns=1000;
mpmat = mpgen(24,3600,1,54321);
usrllh = [40*pi/180 80*pi/180 0];
usrxyz = llh2xyz(usrllh);
loadgps
startt=500; t = startt; deltat=0.1;
randn('state',9083247);
bar1 = waitbar(0,'Calculating Position...   ');
err=[];
for i=1:Ns
    t = t + deltat;
    [svxyzmat,svid,elevation_angle(i,:)] = gensv(usrxyz,t,10,[],6);  % Note the mask angle is set to 2 degrees
    for k=1:4
    sv_xyz{k}(i,:)=( xyz2enu(svxyzmat(k,1:3),usrxyz) )';
    end
    %svxyzmat=svxyzmat(1:end-1,:);
    %svid=svid(1:end-1);
    [prvec,adrvec,err(i,:),common_error(i,:)] = genrng(1,usrxyz,svxyzmat,svid,t,[0.5 0 0 0 0],[],mpmat);
    estusr = olspos(prvec,svxyzmat,usrxyz,10^(-3));  %last argument added by macshen
    
    for k=1:length(svid)  %added by macshen, calculate the true pr error and the residual of the least square fitting
        true_pr_error(i,k)=norm(usrxyz-svxyzmat(k,:))-prvec(k);
        ls_residual(i,k)=norm(estusr(1:3)-svxyzmat(k,:))-prvec(k)+estusr(4);
    end
    true_tilda(i,:)=true_pr_error(i,:)-mean(true_pr_error(i,:));
    
    enuerr(i,:) = ( xyz2enu(estusr(1:3),usrxyz) )';
    terr(i) = estusr(4);  % true clk bias is zero
    waitbar(i/180)
 end
close(bar1)

plot(enuerr(:,1),enuerr(:,2),'*')
axis('equal')
axis('square')
axis([-10 10 -10 10])
grid
title('GPS Positioning Error  -  Static User  -  Half-Hour Scenario')
ylabel('north error (m)')
xlabel('east error (m)')
text(-2,7.5,'Noise, multipath and')
text(-2,6.5,'atmospheric delays simulated')
