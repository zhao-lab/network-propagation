function [ output_args ] = plot_Rand( Position )
%%% Plot the trajectory


TJ_MAX = length(Position);
figure(1);
Margin = 1.1;
ColorOrder = makeColorMap([0 1 0],[1 0 0],[0 0 1],TJ_MAX);


%Plot the trajectory in one map
for i=1:TJ_MAX
    disp(['i = ' num2str(i)]);
    LaLo = [Position{i}.Longitude, Position{i}.Latitude];
    hold on;
    scatter(LaLo(1,1),LaLo(1,2),50,'markeredgecolor','w');
    scatter(LaLo(end,1),LaLo(end,2),50,'d','markeredgecolor','w');
    plot(LaLo(:,1),LaLo(:,2),'linewidth',2,'color',ColorOrder(i,:));
end

%Set the axis limit
xlim = [-83.82,-83.64];ylim = [42.22,42.34];
New_xlim(1) = sum(xlim)/2 - (xlim(2)-xlim(1))/2*Margin;
New_xlim(2) = sum(xlim)/2 + (xlim(2)-xlim(1))/2*Margin;
New_ylim(1) = sum(ylim)/2 - (ylim(2)-ylim(1))/2*Margin;
New_ylim(2) = sum(ylim)/2 + (ylim(2)-ylim(1))/2*Margin;
axis([New_xlim New_ylim]);

xlabel('Longitude(deg)');
ylabel('Latitude(deg)');
set(gca,'fontsize',18);

%Google Map
plot_google_map('MapType','satellite');

end

