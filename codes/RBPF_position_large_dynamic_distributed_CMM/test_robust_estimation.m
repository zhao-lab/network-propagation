%  sndemo05.m       Static User   OLS positioning
clear all
close all
%    
global prvec;
global svxyzmat;  %global quantities added by macshen for testing robust estimation 

mpmat = mpgen(24,3600,1,54321);
usrllh = [0*pi/180 0*pi/180 0];
usrxyz = llh2xyz(usrllh);
loadgps
i=0;
randn('state',9083247);
bar1 = waitbar(0,'Calculating Position...   ');
for t = 41000:1:41200,
    i=i+1;
    [svxyzmat,svid] = gensv(usrxyz,t,2);  % Note the mask angle is set to 2 degrees
    [prvec,adrvec,pr_err_vec(i,:)] = genrng(1,usrxyz,svxyzmat,svid,t,[1 1 0 1 1],[],mpmat);
    estusr = olspos(prvec,svxyzmat);
    fun=@robust_residual;
    robust_estusr=fminunc(fun,estusr);
    enuerr(i,:) = ( xyz2enu(estusr(1:3),usrxyz) )';
    robust_enuerr(i,:) = ( xyz2enu(robust_estusr(1:3),usrxyz) )';
    terr(i) = estusr(4);  % true clk bias is zero
    waitbar(i/180)
end
 norm_err=sqrt(enuerr(:,1).^2+enuerr(:,2).^2);
 robust_norm_err=sqrt(robust_enuerr(:,1).^2+robust_enuerr(:,2).^2);
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



