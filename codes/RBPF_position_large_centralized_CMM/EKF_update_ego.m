function y=EKF_update_ego(pr,svxyzmat,sigma_thermal2,orgxyz);
global EKF;
dx=1;  %here the difference step is chosen large such that floating error will not dominate
sv_n=length(pr);
    pose_local=[EKF.mu(1),EKF.mu(3),0];  %here, the z coordinate is given 0
    pose_xyz=enu2xyz(pose_local,orgxyz);  %give xyz as column vector
    for j=1:sv_n
%         if j<=2
%     range_pre(j,1)=norm(pose_xyz.'-svxyzmat(j,:))+EKF.mu(5)+EKF.mu(6+2*j-1);  %prediction of pr measurement, distance+clock bias+multipath
%         else
    range_pre(j,1)=norm(pose_xyz.'-svxyzmat(j,:))+EKF.mu(5);    
%         end
        end
    
    C=zeros(sv_n,EKF.N);    %generate the observation derivative matrix
    %calculate derivative numerically
    for j=1:sv_n
        pose_xyz_add=enu2xyz(pose_local+[dx,0,0],orgxyz);
        C(j,1)=norm(pose_xyz_add.'-svxyzmat(j,:))-norm(pose_xyz.'-svxyzmat(j,:));
        C(j,1)=C(j,1)/dx;
        pose_xyz_add=enu2xyz(pose_local+[0,dx,0],orgxyz);
        C(j,3)=norm(pose_xyz_add.'-svxyzmat(j,:))-norm(pose_xyz.'-svxyzmat(j,:));
        C(j,3)=C(j,3)/dx; 
        C(j,5)=1;
%         if j<=2
%         C(j,6+2*j-1)=1;
%         end    
    end
    
    innovation=pr.'-range_pre;
    y=individual_compatibility_test(innovation,C,sigma_thermal2);  %check for the compatibility of the measurement, throw away the incompatible measurement
   % y=ones(1,sv_n);
    extract_compatible=[];
    for k=1:length(y)
        if y(k)==1
            extract_compatible=[extract_compatible,k];   %extract the compatible measurement indices
        end
    end
    N_com=length(extract_compatible);
    if N_com~=0   %only update when there is compatible measurement, otherwise, do not update
    C_com=C(extract_compatible,:);
    innovation_com=innovation(extract_compatible);
    Q=eye(N_com,N_com)*sigma_thermal2;
    K=EKF.cov*C_com.'*(Q+C_com*EKF.cov*C_com.')^-1;
    EKF.mu=EKF.mu+K*innovation_com;
    EKF.cov=(eye(EKF.N,EKF.N)-K*C_com)*EKF.cov;
    end
end



function y=individual_compatibility_test(innovation,C,sigma_thermal2)%individual compatibility test
global EKF;
N=length(innovation);
extract=[1,3,5];  %here the clock offset is also accounted for
P_state=EKF.cov(extract,extract);  %extract the submatrix of pose covariance
for k=1:N
    H=C(k,extract);
    P=H*P_state*H.'+sigma_thermal2;
    D2=innovation(k,1)^2/P;
if (D2<=chi2inv(0.98,1))
    y(k)=1;  %indicate compatible
else
    y(k)=0;
end
end  %end for k=1:N
end  %end for function compatibility_test


function y=joint_compatibility_test(innovation,C,sigma_thermal2)
global EKF;
global Best_pair;
global z_JCBB;
global C_JCBB;
z_JCBB=innovation;
C_JCBB=C;
N=length(innovation);
Best_pair=zeros(1,N);
JCBB([],1,sigma_thermal2);
y=Best_pair;
end    


function JCBB(H,i,sigma_thermal2)
global EKF;
global Best_pair;
global z_JCBB;
global C_JCBB;
N=length(z_JCBB);
if (i>N)
    if sum(H>0)>sum(Best_pair>0)
    Best_pair=H;
    end  
else % else for if (i>N)
        if unitary_compatible(i,sigma_thermal2)
        if joint_compatible(H,i,sigma_thermal2)
        JCBB([H 1],i+1,sigma_thermal2);
        end
        end
    if sum(H>0)+N-i>sum(Best_pair>0)
        JCBB([H 0],i+1,sigma_thermal2);
    end
end  %end for if (i>N)
end  %end JCBB


function y=unitary_compatible(i,sigma_thermal2)
global EKF;
global z_JCBB;
global C_JCBB;
extract=[1,3,5];  %here the clock offset is also accounted for
P_state=EKF.cov(extract,extract);  %extract the submatrix of pose covariance
    H=C_JCBB(i,extract);
    P=H*P_state*H.'+sigma_thermal2;
    D2=z_JCBB(i,1)^2/P;
if D2<=chi2inv(0.99,1)
    y=1;
else
    y=0;
end
end


function y=joint_compatible(H,i,sigma_thermal2)
global EKF;
global z_JCBB;
global C_JCBB;
extract=[];
for k=1:length(H)
    if H(k)==1
    extract=[extract,k];
    end
end
extract=[extract,i];  %this is the extracted indices for measurement
extract_pose=[1,3,5];
P_state=EKF.cov(extract_pose,extract_pose);
H=C_JCBB(extract,extract_pose);
CH=H*P_state*H.'+eye(length(extract))*sigma_thermal2;
innovation_JCBB=z_JCBB(extract);
DH2=innovation_JCBB.'*CH^(-1)*innovation_JCBB;       
        if DH2<=chi2inv(0.99,length(extract))  %here N is the total length of H but NN is the # of null hypothesis
            y=1;
        else
            y=0;
        end  %end for if
end  %end for the function