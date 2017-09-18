%  sndemo07.m
%
%     PARITYVEC example:  BIAS error injected onto SV 3 at t = 1750 
%                        (750 seconds after start)
%
clear all
close all
%    
mpmat=mpgen(24,3600);
usrllh = [0*pi/180 0*pi/180 0];
usrxyz = llh2xyz(usrllh);
loadgps
i=0;  initpos = [0 0 0 0];
randn('state',309874);
bar1 = waitbar(0,'Finding Residuals...  ');
for t = 1000:5:1750,  %changed by macshen
    i=i+1;
    [svxyzmat,svid] = gensv(usrxyz,t);
    [prvec,adrvec] = genrng(1,usrxyz,svxyzmat,svid,t,[1 1 0 1 1],[],mpmat);
    if t > 1750,
       prvec(1) = prvec(1) + 200;
    end
    estusr = olspos(prvec,svxyzmat);
    enuerr(i,:) = ( xyz2enu(estusr(1:3),usrxyz) )';  %added by macshen
    sv(i) = svid(1);
    parvec = parityvec(prvec,svxyzmat,estusr);
    p(i) = norm(parvec);
    time(i) = t;
    terr(i) = estusr(4);  % true clk bias is zero
    waitbar(i/360)
 end
 close(bar1);
plot(time,p)
title('Fault Detection Example: Bias on SV 3 Injected at t = 1750')
ylabel('magnitude of parity vector')
xlabel('GPS time of week in seconds')
