% bacteria post processor
Run_Name = 'Figure1_Combined';
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


for k = 1:6
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

% Set up MSD Vector for Each Time Step
[nTimeSteps, nCells, nSimulations] = size(b_trajectoriesX);
MSD = zeros(nTimeSteps,nSimulations);

for i = 1:nSimulations
    MSD_Simulation = zeros(nTimeSteps,nCells);
    for j=1:nCells
        MSD_Simulation(:,j) = (b_trajectoriesX(:,j,i) - b_trajectoriesX(1,j,i)).^2 + (b_trajectoriesY(:,j,i) - b_trajectoriesY(1,j,i)).^2;
    end
    MSD(:,i) = mean(MSD_Simulation,2);
end

Diffusion = MSD;
for t = 2:nTimeSteps
    Diffusion(t,:) = Diffusion(t,:)/(4*.1*t);
end

sourcesX = 2000; %C{1} ;
sourcesY = 2000; %C{2} ;

distToSourceX = [b_trajectoriesX] - sourcesX ;
distToSourceY = [b_trajectoriesY] - sourcesY ;
distToSource = ( distToSourceX .* distToSourceX + distToSourceY .* distToSourceY ).^ 0.5 ;
xDrift2 = distToSource(1,:,:) - [distToSource] ;
vDrift2 = xDrift2(2:end,:,:) ./  ([2:2001]'./10);
% 
vdAvg = [mean(vDrift2,2), std(vDrift2,0,2)] ;
distToSourceAvg = [mean(distToSource,2), std(distToSource,0,2)] ;

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

save('Figure1_Workspace')

%% 1 Pannel 
t = tiledlayout(1,1);
t.TileSpacing = 'compact'; % Reduces space between tiles
t.Padding = 'compact';     % Reduces space around the outside of the tiles

% Initialize an array to store plot handles for the legend
plotHandles = gobjects(6, 1);

%% Plot Average Distance to Source
ax1 = nexttile;
line_color = {'#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b'};
for i = 1:6
    plotHandles(i) = plot(1/10:.1:2001/10, distToSourceAvg(:,1,i), 'Color', line_color{i},'LineWidth',2);
    hold on
end  
xlim([1 200]);
xlims = xlim;
% Define the y-value for the dashed line
y_value = distToSourceAvg(1,1,1);

% Add a dashed line at y = y_value
line(xlims, [y_value y_value], 'LineStyle', '--', 'Color', 'k','LineWidth',2); % '--' makes the line dashed, 'Color', 'r' makes the line red. Adjust as needed.
xlabel('Time (s)','FontSize',14);
ylabel('Distance from the Center (\mum)','FontSize',14);
ax1.FontSize = 14;
%text(ax1, -0.1, 1.00, 'A', 'Units', 'normalized', 'FontSize', 18, 'FontWeight', 'bold');

% Create a single legend for all plots outside the tiled layout
l = legend(plotHandles, {'No Wrap', ['Wrap = 15' char(176)], ['Wrap = 45' char(176)], ['Wrap = 90' char(176)], ['Wrap = 135' char(176)], ['Wrap = 180' char(176)]},'Location', 'northwest');
% l = legend(plotHandles, {'No Wrap', 'Wrap = 15\u00B0', 'Wrap = 45\u00B0', 'Wrap = 90\u00B0', 'Wrap = 135\u00B0', 'Wrap = 180\u00B0'});
l.FontSize = 14;

% Adjust the figure size if necessary
set(gcf,'units','inches','pos',[0 0 7 6])
pos = get(gcf, 'Position');
set(gcf, 'PaperPositionMode', 'Manual', 'PaperUnits', 'inches', 'PaperSize', [pos(3) pos(4)], 'PaperPosition', [0 0 pos(3) pos(4)]);

%title(t,'Bacterial Migration With Varying Wrap Probabilities', 'FontSize', 14, 'FontWeight', 'bold')

% % Save the figure
dpi = 600;

% Define the base filename (without extension)
baseFilename = fullfile(pwd, Run_Name, '1_6_Average_Distance');

% Save as .fig
savefig([baseFilename, '.fig']);

% Save as .pdf (with specified DPI)
print(gcf, [baseFilename, '.pdf'], '-dpdf', ['-r', num2str(dpi)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all

% Max_Drift = max(vdAvg,[],1);
% for j = 1:6 
%     disp(['Simulation' num2str(j) ])
%     disp(['Max Drift Velocity is ' num2str(Max_Drift(1,1,j)) ])
%     disp(['Ending Aggregation Time is ' num2str(tIn_avg(end,1,j)) ])
% end

for i = 1:6
    for k = 1:5
        figure(k)
        % Plot the trajectory in red
        plot([b_trajectoriesX(:,k,i)-2000], [b_trajectoriesY(:,k,i)-2000], 'b-', LineWidth=1)
        hold on        
        % Add a black dot at the start of the trajectory
        plot(b_trajectoriesX(1,k,i)-2000, b_trajectoriesY(1,k,i)-2000, 'ko', 'MarkerFaceColor', 'r', 'MarkerSize', 4)
        % Add a blue dot at the end of the trajectory
        plot(b_trajectoriesX(end,k,i)-2000, b_trajectoriesY(end,k,i)-2000, 'bo', 'MarkerFaceColor', 'k', 'MarkerSize', 4)
        hold off
        % Label both x and y axes as "Microns"
        xlabel('Microns (\mum)','FontSize',14);
        ylabel('Microns (\mum)','FontSize',14);
        set(gca, 'FontSize', 14);
        set(gcf,'units','inches','pos',[0 0 7 6])
        pos = get(gcf, 'Position');
        set(gcf, 'PaperPositionMode', 'Manual', 'PaperUnits', 'inches', 'PaperSize', [pos(3) pos(4)], 'PaperPosition', [0 0 pos(3) pos(4)]);
        % % Save the figure
        dpi = 600;
        % Define the base filename (without extension)
        baseFilename = fullfile(pwd, strcat('/', Run_Name,'/Trajectory', num2str(i), '-', num2str(k)));
        % Save as .fig
        savefig([baseFilename, '.fig']);
        % Save as .pdf (with specified DPI)
        print(gcf, [baseFilename, '.pdf'], '-dpdf', ['-r', num2str(dpi)]);
    end
end