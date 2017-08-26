clearvars -except vehicle

R=3000;
con=zeros(50);
for k=1:length(vehicle(:,1))
    for j=1:length(vehicle(:,1))
        r=norm(vehicle(k,1:2)-vehicle(j,1:2));
        if r<R
            con(k,j)=1;
        end
    end
end

s_con=sum(con');
de_node=find(s_con>=18);

m=0;
for k=1:50
    if ~ismember(k,de_node)
        new_vehicle(m+1,1:3)=vehicle(k,1:3);
        m=m+1;
    end
end

for k=1:length(new_vehicle(:,1))
    for j=1:length(new_vehicle(:,1))
        r=norm(new_vehicle(k,1:2)-new_vehicle(j,1:2));
        if r<R
            new_con(k,j)=1;
        end
    end
end

a=sum(new_con');
histogram(a)
