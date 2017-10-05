clearvars -except vehicle
load medium_sparse.mat

% Added for test:
vehicle = new_vehicle;
R=3000;
con=zeros(size(new_vehicle,1));
for k=1:length(vehicle(:,1))
    for j=1:length(vehicle(:,1))
        r=norm(vehicle(k,1:2)-vehicle(j,1:2));
        if r<R
            con(k,j)=1;
        end
    end
end

s_con=sum(con');
% delete any nodes that not within the distance range of 18?
de_node=find(s_con>=14);

m=0;
for k=1:size(new_vehicle,1)
    if ~ismember(k,de_node)
        new_vehicle1(m+1,1:3)=vehicle(k,1:3);
        m=m+1;
    end
end

for k=1:length(new_vehicle1(:,1))
    for j=1:length(new_vehicle1(:,1))
        r=norm(new_vehicle1(k,1:2)-new_vehicle1(j,1:2));
        if r<R
            new_con(k,j)=1;
        end
    end
end

a=sum(new_con');
figure;
histogram(a)
% Is that true?
xlabel('Value of cons distance')
ylabel('Number of value appearance')
