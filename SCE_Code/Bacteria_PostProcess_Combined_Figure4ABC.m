% bacteria post processor
Run_Name = 'Figure4ABC_Combined';
mkdir (Run_Name);

nbacteria = 50 ;
domainSize = 1000 ;

b_trajectoriesX = [] ;
b_trajectoriesY = [] ;
b_locVelocity = [] ;
b_locFriction = [] ;
b_Orientation = [] ;
b_chemical = [] ;

sourcesX = [] ;
sourcesY = [] ;

nReverse = [] ;
time_toSource = [] ;
time_inSource = [] ;
init_maxRunDuration = [] ;

vDrift = [];

R_threshold = 50 ;
t_in = [] ;
t_to = [] ;


for k = 1:2
    for fileNum = 1:2001
        % disp(fileNum)

        filename = [pwd,'/dataStats/machine',num2str(k), '/Bacteria_' num2str(fileNum) '.txt'];
        fileID = fopen(filename);
        C = textscan(fileID,'%f %f %f %f %f %f %f','Delimiter','\t');
        fclose(fileID);
        b_trajectoriesX(fileNum,:,k) = [transpose( C{1}) ] ;
        b_trajectoriesY(fileNum,:,k) = [transpose( C{2}) ] ;
        b_locVelocity(fileNum,:,k)   = [transpose( C{3}) ] ;
        b_locFriction(fileNum,:,k)   = [transpose( C{4}) ] ;
        b_Orientation(fileNum,:,k)   = [transpose( C{5}) ] ;
        b_chemical(fileNum,:,k)      = [transpose( C{6}) ] ;

    end
end
sourcesX = 2000; %C{1} ;
sourcesY = 2000; %C{2} ;

distToSourceX = [b_trajectoriesX] - sourcesX ;
distToSourceY = [b_trajectoriesY] - sourcesY ;
distToSource = ( distToSourceX .* distToSourceX + distToSourceY .* distToSourceY ).^ 0.5 ;

Displacement = sqrt((b_trajectoriesX(2:end,:,:) - b_trajectoriesX(1:end-1,:,:)).^2 + (b_trajectoriesY(2:end,:,:) - b_trajectoriesY(1:end-1,:,:)).^2);
Cumulative_Displacement = cumsum(Displacement,1);

%% Drift Velocity by Distance
xDrift2 = distToSource(1,:,:) - [distToSource] ;
vDrift2 = xDrift2(2:end,:,:) ./  (Cumulative_Displacement);

vdAvg = [mean(vDrift2,2), std(vDrift2,0,2)] ;
distToSourceAvg = [mean(distToSource,2), std(distToSource,0,2)] ;
distToSource_Lower_STD(:,1,1) = distToSourceAvg(:,1,1)-distToSourceAvg(:,2,1);
distToSource_Upper_STD(:,1,1) = distToSourceAvg(:,1,1)+distToSourceAvg(:,2,1);
distToSource_Lower_STD(:,1,2) = distToSourceAvg(:,1,2)-distToSourceAvg(:,2,2);
distToSource_Upper_STD(:,1,2) = distToSourceAvg(:,1,2)+distToSourceAvg(:,2,2);

dist2D = [sqrt((b_trajectoriesX - sourcesX).^2 + (b_trajectoriesY - sourcesY).^2 ) ] ;
t_in = zeros(size(dist2D,1),size(dist2D,2),size(dist2D,3));

for l = 1:size(dist2D,3)
    for i = 1:size(dist2D,1)
        for j = 1:size(dist2D,2)

            if (dist2D(i,j,l)< R_threshold )

                t_in(i+1,j,l) = t_in(i,j,l)+ 0.1 ; 
            else
                t_in(i+1,j,l) = t_in(i,j,l) ;

            end

        end
    end
end

tIn_avg = [mean(t_in,2)] ;
save('Figure4ABC_Workspace')

%% Create a tiled layout
t = tiledlayout(1,3);
t.TileSpacing = 'compact'; % Reduces space between tiles
t.Padding = 'compact';     % Reduces space around the outside of the tiles

% Initialize an array to store plot handles for the legend
plotHandles = gobjects(2, 1);

%% Plot Average Distance to Source
ax1 = nexttile;
line_color = {'#1f77b4','#2ca02c'};
xcoords = 1/10:.1:2001/10;
i = 1;
    fill([xcoords, flip(xcoords)], [distToSource_Upper_STD(:,1,i); flip(distToSource_Lower_STD(:,1,i))], [0.7, 0.85, 1])
    hold on
i = 2;
    fill([xcoords, flip(xcoords)], [distToSource_Upper_STD(:,1,i); flip(distToSource_Lower_STD(:,1,i))], [0.7, 1, 0.7])
    hold on

for i = 1:2
    plotHandles(i) = plot(1/10:.1:2001/10, distToSourceAvg(:,1,i), 'Color', line_color{i},'LineWidth',3);
    hold on
end  

xlim([1 200]);
xlims = xlim;
% Define the y-value for the dashed line
y_value = distToSourceAvg(1,1,1);

% Add a dashed line at y = y_value
line(xlims, [y_value y_value], 'LineStyle', '--', 'Color', 'k','LineWidth',2); % '--' makes the line dashed, 'Color', 'r' makes the line red. Adjust as needed.
ax1.FontSize = 20;
xlabel('Time (s)','FontSize',28);
ylabel('Avg. Dist. to Source (\mum)','FontSize', 28);
text(ax1, -0.17, 1.00, 'A', 'Units', 'normalized', 'FontSize', 28, 'FontWeight', 'bold');

%% Plot Drift Velocity
ax2 = nexttile;
for i = 1:2
    plot(1/10:.1:2000/10, vdAvg(:,1,i), 'Color', line_color{i},'LineWidth',3);
    hold on
end 
xlim([1 200])
xlims = xlim;
% Define the y-value for the dashed line
y_value = 0;

% Add a dashed line at y = y_value
line(xlims, [y_value y_value], 'LineStyle', '--', 'Color', 'k','LineWidth',3); % '--' makes the line dashed, 'Color', 'r' makes the line red. Adjust as needed.
ax2.FontSize = 20;
xlabel('Time (s)','FontSize',28) 
ylabel('Directional Efficiency (%)','FontSize',28)
text(ax2, -0.13, 1.00, 'B', 'Units', 'normalized', 'FontSize', 28, 'FontWeight', 'bold');


%% Plot Aggregation Time
ax3 = nexttile;
for i = 1:2
    plot(1/10:.1:2002/10, tIn_avg(:,1,i), 'Color', line_color{i},'LineWidth',3);
    hold on
end 
xlim([1 200])
ax3.FontSize = 20;
xlabel('Time (s)','FontSize',28) 
ylabel('Aggregation Time (s)','FontSize',28)
text(ax3, -0.15, 1.00, 'C', 'Units', 'normalized', 'FontSize', 28, 'FontWeight', 'bold');

% Create a single legend for all plots outside the tiled layout
l = legend(ax3, plotHandles, {'Prob = 0%', 'Prob = 38%'});
l.Location = 'northwest';
l.FontSize = 26;
% Adjust the figure size if necessary
set(gcf,'units','inches','pos',[0 0 24 7])
set(gcf, 'Units', 'inches');
pos = get(gcf, 'Position');
set(gcf, 'PaperPositionMode', 'Manual', 'PaperUnits', 'inches', 'PaperSize', [pos(3) pos(4)], 'PaperPosition', [0 0 pos(3) pos(4)]);

% Construct the file path for PDF
filename = fullfile(pwd, Run_Name, 'Combined_Plots.pdf');

% Save as PDF
print(gcf, filename, '-dpdf', '-vector','-r600');


for i = 1:2
    for j = 1:800
        slope(i,j) = (distToSourceAvg(j+1,1,i)-distToSourceAvg(j,1,i))/.1;
    end
end
avg_slope = mean(slope,2);

%% Average Distance Difference 
dist_difference = distToSourceAvg(:,1,1) - distToSourceAvg(:,1,2);
max_difference = max(dist_difference);

%% VDrift Difference 
Vdrift_difference = vdAvg(:,1,2) - vdAvg(:,1,1);
max_difference_Vdrift = max(Vdrift_difference);

% 
% Max_Drift = max(vdAvg,[],1);
% for j = 1:6 
%     disp(['Simulation' num2str(j) ])
%     disp(['Max Drift Velocity is ' num2str(Max_Drift(1,1,j)) ])
%     disp(['Ending Aggregation Time is ' num2str(tIn_avg(end,1,j)) ])
% end

% for k = 1:5
%     figure(k)
%     plot([b_trajectoriesX(:,k,7)-500],[b_trajectoriesY(:,k,7)-500])
%     saveas(gcf,[pwd, strcat('/', Run_Name,'/MWCTrajectory') num2str(k) '.png'])
% end
