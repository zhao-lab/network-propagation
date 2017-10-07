function [svxyzmat,svid,elevation_angle,sigma2] = gensv(usrxyz,time,maskang,orgxyz,sv_number)
%GENSV	Generate matrix of positions of satellites which are
%	visible to the user at location usrxyz.  If orgxyz is not
%	specified, then satellite positions are computed in ECEF 
%   coordinates.  If orgxyz is specified, it is taken to be 
%   the origin of a locally-level, East-North-UP (ENU) 
%   coordinate system and satellite positions are given
%	accordingly.
%
%	svxyzmat = GENSV(usrxyz,time,maskang,orgxyz)
%                       or
%       [svxyzmat,svid] = GENSV(usrxyz,time,maskang,orgxyz)
%
%   INPUTS
%	usrxyz(1:3) = user position in cartesian ECEF coordinates if
%	              orgxyz is not provided.  If
%	              orgxyz is provided, usrxyz is given in cartesian
%	              ENU coordinates relative to orgxyz.
%	time = time of week in seconds
%	maskang = satellite visibility mask angle (in degrees).  This
%	          argument is optional and has a default value of 5 degrees
%	orgxyz(1:3) = optional argument.  Cartesian ECEF coordinates of
%	              locally-level ENU coordinate origin
%   sv_number specify the wanted number of satellites, added by macshen

%   OUTPUTS
%	svxyzmat = matrix of visible satellite positions in cartesian 
%                  coordinates.  
%                  svxyzmat(i,1:3) = x,y,z coordinates for satellite i
%   svid = vector of identification numbers for the visible satellites
%              corresponding to svxyzmat
%
%   NOTE: GENSV expects satellite orbital parameters to be loaded and maintained
%         as global variables.  This is accomplished with LOADGPS, LOADGLO or LOADGG.

%	M. & S. Braasch 11-96
%	Copyright (c) 1996 by GPSoft
%	All Rights Reserved.
%

global SVIDV MV OMGV RV INCLV TOEV

if nargin<4,orgxyz=[];end
if nargin<3,maskang=5;end
if nargin<2,error('insufficient number of input arguments');end
if isempty(orgxyz),usrece=usrxyz;else,usrece=enu2xyz(usrxyz,orgxyz);end
if isempty(orgxyz),enuflg=0;else,enuflg=1;end
%%
i = 0;
numsv = max(size(SVIDV));
svid = 0;
for N = 1:numsv
   svxyz = svpos(RV(N),TOEV(N),MV(N),OMGV(N),INCLV(N),time); 
   svenu = xyz2enu(svxyz,usrece);  % User is origin
   el = (180/pi)*atan2(svenu(3),norm(svenu(1:2)));   
   if el > maskang
      i = i + 1;
      if enuflg==0
	svxyzmat(i,:) = svxyz';      % ECEF position
      else
	svxyzmat(i,:) = (xyz2enu(svxyz,orgxyz))';      % ENU position
      end
      svid(i) = SVIDV(N);
      elevation_angle(i)=el;%added by macshen, to output the elevation angle
      sigma2(i)=16+9/sin(elevation_angle(i)/180*pi)^2; %added by macshen, 
      %variance of the noise, reference Single-point Positioning with the Pseudorange of Single-frequency GPS Considering the Stochastic Model
   end
end

%%

if nargin<5
    sv_number=length(svid);
else  %noted that there can be case where the wanted number is larger than the available number of satellites
    svid=svid(1:sv_number);
    svxyzmat=svxyzmat(1:sv_number,:);
    elevation_angle=elevation_angle(1:sv_number);
    sigma2=sigma2(1:sv_number);
end
