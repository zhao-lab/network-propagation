function [estusr,H] = olspos_clone(prvec,svxyzmat,v_id,orgxyz,initpos,tol,weight_variance)  %modified by macshen, now the matrix H is also output, and weight has been used for WLS
%OLSPOS	Compute position from satellite positions and pseudoranges
%       via ordinary least squares.
%
%	estusr = OLSPOS(prvec,svxyzmat,initpos,tol)
%
%   INPUTS
%	prvec = vector of 'measured' pseudoranges for satellites
%               specified in svxyzmat
%	svxyzmat(i,1:3) = position of satellite i in user defined 
%                     cartesian coordinates.
%	initpos = optional argument.  Initial 'estimate' of user state:
%                 three-dimensional position and clock offset
%                 (in user defined coordinates).  Used to speed up 
%                 iterative solution.  Initial clock offset is
%                 optional with default value = 0.
%	tol = optional argument.  Tolerance value used to determine 
%             convergence of iterative solution.  Default value = 1e-3
%
%   OUTPUTS
%	estusr(1:3) = estimated user x, y, and z coordinates
%	estusr(4) = estimated user clock offset
%          Note: all four elements of estusr are in the same units
%                as those used in prvec

%	M. & S. Braasch 11-96
%	Copyright (c) 1996 by GPSoft
%	All Rights Reserved.
%
nn=length(prvec);
% if v_id==3
%     prvec=prvec(1:nn-1);
%     svxyzmat=svxyzmat(1:nn-1,:);
% end
    
if nargin<6,tol=1e-3;end
if nargin<5,initpos=[0 0 0 0];end
if nargin<2,error('insufficient number of input arguments'),end
[m,n]=size(initpos);
if m>n, estusr=initpos';else,estusr=initpos;end
if max(size(estusr))<3
   error('must define at least 3 dimensions in INITPOS')
end
if max(size(estusr))<4,estusr=[estusr 0];end
numvis=max(size(svxyzmat));
beta=[1e9 1e9 1e9 1e9];
maxiter=10;
iter=0;
while ((iter<maxiter)&&(norm(beta)>tol))
    for N = 1:numvis
	pr0 = norm(svxyzmat(N,:)-estusr(1:3));
	y(N,1) = prvec(N) - pr0 - estusr(4);
    end
    
    H = hmat(svxyzmat,estusr(1:3)); %hamt: outputs direction cosine matrix for GPS positioning
    % check "help hmat"     
    %the direction cosines (or directional cosines) of a vector are the cosines of the angles between the vector and the three coordinate axes. 

    if nargin<5
       beta = H\y;
    else
    R_minus1=diag(1./weight_variance);
    A=H.'*R_minus1*H;
    beta=A\(H.'*R_minus1*y);     %the former three lines are added by macshen for WLS
    end
    
    estusr=estusr+beta';
    iter=iter+1;
end

for k=1:length(svxyzmat(:,1))
    svxyzmat(k,1:3)=xyz2enu(svxyzmat(k,1:3),orgxyz);
end
H=hmat(svxyzmat,[0,0,0]);
