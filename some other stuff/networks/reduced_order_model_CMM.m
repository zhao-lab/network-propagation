%this script simulates a heuristic reduced order model for the CMM of four vehicle 
clear all
x=zeros(8,1);   %each vehicle has two scalar as x and y error
x(1,1)=20;
x(3,1)=-100
w_direction=[1;0;-1;0;0;1;0;-1];   %the direction of Gaussian random driving force
n_step=25;   %number of fusion step
w_mag=sqrt(n_step)*0.1;  %magnitude of the driving force
for k=1:300
    x(:,k+1)=x(:,k)+0*w_mag*abs(randn(8,1)).*w_direction;
    
    s12_input=[-(x(1,k+1)-abs(x(1,k+1)))^2/4-x(2,k+1)^2;-(x(3,k+1)-abs(x(3,k+1)))^2/4-x(4,k+1)^2];
    s23_input=[-(x(3,k+1)+abs(x(3,k+1)))^2/4-x(4,k+1)^2;-(x(5,k+1)+abs(x(5,k+1)))^2/4-x(6,k+1)^2];
    s34_input=[-(x(6,k+1)-abs(x(6,k+1)))^2/4-x(5,k+1)^2;-(x(8,k+1)-abs(x(8,k+1)))^2/4-x(7,k+1)^2];
    s41_input=[-(x(8,k+1)+abs(x(8,k+1)))^2/4-x(7,k+1)^2;-(x(2,k+1)+abs(x(2,k+1)))^2/4-x(1,k+1)^2];
        s12=softmax(s12_input);
        s23=softmax(s23_input);
        s34=softmax(s34_input);
        s41=softmax(s41_input);
        A=[s12(1),0,s12(2),0,0,0,0,0;
           0,s12(1),0,s12(2),0,0,0,0;
           0,0,s23(1),0,s23(2),0,0,0;
           0,0,0,s23(1),0,s23(2),0,0;
           0,0,0,0,s34(1),0,s34(2),0;
           0,0,0,0,0,s34(1),0,s34(2);
           s41(2),0,0,0,0,0,s41(1),0;
           0,s41(2),0,0,0,0,0,s41(1);];
       x(:,k+1)=A*x(:,k+1);
end

           
        