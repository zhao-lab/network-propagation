function enu = xyz2enu_exp(xyz)%,orgxyz) only need one input%by macshen
%XYZ2ENU	Convert from WGS-84 ECEF cartesian coordinates to 
%               rectangular local-level-tangent ('East'-'North'-Up)
%               coordinates.
%
%	enu = XYZ2ENU(xyz,orgxyz)	
%
%    INPUTS
%	xyz(1) = ECEF x-coordinate in meters
%	xyz(2) = ECEF y-coordinate in meters
%	xyz(3) = ECEF z-coordinate in meters
%
%	orgxyz(1) = ECEF x-coordinate of local origin in meters
%	orgxyz(2) = ECEF y-coordinate of local origin in meters
%	orgxyz(3) = ECEF z-coordinate of local origin in meters
%
%    OUTPUTS
%       enu:  Column vector
%		enu(1,1) = 'East'-coordinate relative to local origin (meters)
%		enu(2,1) = 'North'-coordinate relative to local origin (meters)
%		enu(3,1) = Up-coordinate relative to local origin (meters)

%	Reference: Alfred Leick, GPS Satellite Surveying, 2nd ed.,
%	           Wiley-Interscience, John Wiley & Sons, 
%	           New York, 1995.
%
%	M. & S. Braasch 10-96
%	Copyright (c) 1996 by GPSoft
%	All Rights Reserved.

% if nargin<2,error('insufficient number of input arguments'),end
%tmpxyz = xyz;
% tmporg = orgxyz;
% if size(tmpxyz) ~= size(tmporg), tmporg=tmporg'; end,
% difxyz = tmpxyz - tmporg;
% [m,n] = size(difxyz); if m<n, difxyz=difxyz'; end,
% orgllh = xyz2llh(orgxyz);
% phi = orgllh(1);
% lam = orgllh(2);
% sinphi = sin(phi);
% cosphi = cos(phi);
% sinlam = sin(lam);
% coslam = cos(lam);
% R = [ -sinlam          coslam         0     ; ...
%       -sinphi*coslam  -sinphi*sinlam  cosphi; ...
%        cosphi*coslam   cosphi*sinlam  sinphi];
% enu = R*difxyz;
ref_llh=[42.29724065/180*pi -83.7065695/180*pi 265.0631];
llh=xyz2llh(xyz);
dif_llh=llh-ref_llh;
x=dif_llh(2)*6371000*cos(ref_llh(1));
y=dif_llh(1)*6371000;
z=dif_llh(3)-0.0273*y;  %here xyz is the local coordinate in the plane

if abs(x)>5000|abs(y)>5000|abs(z)>2500   %not local transformation, large error expected
    'error line 56 xyz2enu_exp, global transformation'
    a  %force the program to stop
end
enu=[x;y;z];
end

