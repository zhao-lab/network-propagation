function [x_star,y_star] = prestige(A)
%routine by Mo-Han Hseih to calculate prestige and acquaintance based on
%adjacency matrix A
%usage [x,y]=prestige(A); x is prestige, y is acquaintance
AAT=A * A';
ATA=A' * A;
[V1 D1]=eig(ATA);
m1=size(V1,2);
x_star=V1(:,m1);
[V2 D2]=eig(AAT);
m2=size(V2,2);
y_star=V2(:,m2);