%  sndemo10.m
%
%  Demonstrate vehicle path generator
%
clear all
close all
%    
mpmat=mpgen(24,3600,1,54321);
orgllh = [0*pi/180 0*pi/180 0];
orgxyz = llh2xyz(orgllh);
loadgps
startt=1000; t = startt; deltat=5;
segp = [150 90 .2; 150 90 .2; 150 90 .2];
usrenu = pathgen([0 0 0],[5 0],segp,deltat);
EndLoop = max(size(usrenu));
bar1 = waitbar(0,'Generating User Solutions...  ');
%id=[]; %added by macshen
%pr_err=[]; %added by macshen
for i = 1:EndLoop,
    t = t + deltat;
    usrxyz=enu2xyz(usrenu(i,:),orgxyz);
    [svxyzmat,svid] = gensv(usrxyz,t,5);
   % id=[id;svid(1)];  %added by macshen
    [prvec,adrvec,pr_err_vec] = genrng(1,usrxyz,svxyzmat,svid,t,[1 1 0 1 1],[],mpmat);  %pr_err_vec added by macshen
   % pr_err=[pr_err;pr_err_vec(1)]; %added by macshen
    estusr = olspos(prvec,svxyzmat);
    estenu(i,:) = ( xyz2enu(estusr(1:3),orgxyz) )';
    err(i,1:3) = estenu(i,1:3) - usrenu(i,:);
    terr(i) = estusr(4);  % true clk bias is zero
 	 waitbar(i/EndLoop)   
  end
close(bar1);
plot(usrenu(:,1),usrenu(:,2),'-',estenu(:,1),estenu(:,2),'*')
title('True and Estimated Trajectories')
ylabel('north direction [meters]')
xlabel('east direction [meters]')

