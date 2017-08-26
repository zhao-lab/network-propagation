clear all
for tt=1:10
clearvars -except tt a er e_max
d=0.75;
mag=0.5;
N=500;
x=zeros(4,1);
y=zeros(4,1);
    com_list{1}=[1,2];
    com_list{2}=[2,3];
    com_list{3}=[3,4];
    com_list{4}=[4,1];
for k=1:N
    dx=rand(4,1)*mag-mag/2;
    dy=rand(4,1)*mag-mag/2;
    if x(1,k)<-d
        dx(1)=dx(1)+mag/2;
    end
    if y(1,k)<-d
        dy(1)=dy(1)+mag/2;
    end
    if x(2,k)>d
        dx(2)=dx(2)-mag/2;
    end
    if y(2,k)<-d
        dy(2)=dy(2)+mag/2;
    end
    if x(3,k)>d
        dx(3)=dx(3)-mag/2;
    end
    if y(3,k)>d
        dy(3)=dy(3)-mag/2;
    end
    if x(4,k)<-d
        dx(4)=dx(4)+mag/2;
    end
    if y(4,k)>d
        dy(4)=dy(4)-mag/2;
    end
    dx_his(:,k)=dx;
    dy_his(:,k)=dy;
    x(:,k+1)=x(:,k)+dx;
    y(:,k+1)=y(:,k)+dy;
    
% 
%     com_err=[x(1,k+1),x(2,k+1),x(3,k+1),x(4,k+1);y(1,k+1),y(2,k+1),y(3,k+1),y(4,k+1)];
%     mean_com=mean(com_err')';
%     for j=1:4
%     com_err(1:2,j)=com_err(1:2,j)-mean_com;
%     end
% for k1=1:4
% for k2=1:4
% M_opt(k1,k2)=com_err(:,k1).'*com_err(:,k2);  %matrix fed to the optimization program
% end
% end
% 
% for j=1:4
% N_alloc{j}=opt_com_weight(M_opt,j,50,com_list{j});   %calculate the number of particles for allocation
% N_alloc{j}=N_alloc{j}/sum(N_alloc{j});
% end
% 
% A=[N_alloc{1}';N_alloc{2}';N_alloc{3}';N_alloc{4}'];
    
alpha=0.95;
    A=[alpha,1-alpha,0,0;
        0,alpha,1-alpha,0;
        0,0,alpha,1-alpha;
        alpha,0,0,1-alpha];
    
    
    x(:,k+1)=A*x(:,k+1);
    y(:,k+1)=A*y(:,k+1);
end

    a(tt)=sqrt(sum((y(1,:)-y(2,:)).^2+(x(1,:)-x(2,:)).^2)/N)
    er(tt)=sqrt((sum(x(1,:).^2)+sum(y(1,:).^2))/N)
    e_max(tt)=sqrt(max(x(1,:).^2+y(1,:).^2));
    end
    
    mean(a)
    mean(er)
    mean(e_max)
    