function w=robust_residual(x)  
%x is the vector containing position and distance caused by delay time
%created by macshen
global prvec;
global svxyzmat;

N=length(prvec);
for k=1:N
    w(k)=(x(1)-svxyzmat(k,1))^2+(x(2)-svxyzmat(k,2))^2+(x(3)-svxyzmat(k,3))^2;
    w(k)=sqrt(w(k))-x(4)-prvec(k);
    w(k)=sqrt(w(k)^2/2+1)-1;
    %w(k)=abs(w(k));
end
w=sum(w);
end