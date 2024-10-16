
% load data
clear
load('data.mat')
%trace{1}.sign_winding

% Select traces with "good" count traces being selected
count = 0;
for i = 1:numel(trace)
    if strcmp(trace{i}.comments,'good')
    count = count + 1;
    end
end

% Extract parameters (height, force, topo speed, and wait time, and name of these traces) from the selected traces with "good" comment
time_holding_1 = [];
turn_holding_1 = [];
z_holding_1 = [];
time_holding_2 = [];
z_holding_2 = [];
height = zeros(count, 1);
force = zeros(count,1);
%names = {};


index = 1;
for i = 1 : numel(trace)
    if strcmp(trace{i}.comments,'good')
       time_holding_1(:,index) = trace{i}.time_holding1;
       turn_holding_1(:,index) = trace{i}.turn_holding1;
       time_holding_2(:,index) = trace{i}.time_holding2;
       temp_height_0 = mean(trace{i}.z_holding1(trace{i}.time_holding1<=10));
       z_holding_1(:,index) = trace{i}.z_holding1/temp_height_0;%normolized height
       z_holding_2(:,index) = trace{i}.z_holding2/temp_height_0;%normolized height
       height(index)= temp_height_0;
       force(index) = trace{i}.F_from_var_length_p;
       index = index +1;
    end   
end
index_height = find(height>2.9);
height = height(index_height);
time_holding_1 = time_holding_1(:, index_height);
%time_holding_2 = time_holding_2(:, index_height);
z_holding_1 = z_holding_1(:, index_height);
z_holding_2 = z_holding_2(:, index_height);
turn_holding_1 = turn_holding_1(:,1);

Exten_end = zeros(numel(index_height),0);
% Create figure
figure1 = figure;
% Create axes
axes1 = axes('Parent',figure1,...
    'Position',[0.13 0.11 0.542916666666667 0.815]);
hold(axes1,'on');
for i = 1:numel(index_height)
    time_holding = time_holding_1(:,1);
    t_temp = time_holding;
%     rate = 1/0.4;  
      dt = mean(diff(time_holding));
%     window = round(1/rate / dt);
     order = 2;
%     if mod(window, 2) == 0
%            window = window+1;
%     end
%     Z_holding_s1 = sgolayfilt(z_holding_1(:,i),order,window);
         
    window2 = round(5/ dt);
    if mod(window2, 2) == 0
           window2 = window2+1;
    end
    Z_holding_s2 = sgolayfilt(z_holding_1(:,i),order,window2);
    Exten_end(i) = mean(Z_holding_s2(end-0.1*100:end));
    plot1 = plot(t_temp, Z_holding_s2,'Color',[0.800000011920929 0.800000011920929 0.800000011920929]);
    %plot(time_holding_1(:,1),z_holding_1(:,i))
end
Z_holding_avg = nanmean(z_holding_1,2);
t = time_holding_1(:,1);
plot(t, Z_holding_avg);
%turn_index = find(turn_holding_1 == 50);
%t_0 = t(turn_index(1))+10;
%line([t(turn_index(1)) t(turn_index(1))], [1.5 4],'LineStyle','--','LineWidth',1, 'Color','k');
line([ 0 1000], [1 1],'LineStyle','--','LineWidth',1, 'Color','k');
t_0 = 90;
% t_index = t(turn_index(1)).
%hold on;
Z_holding_avg_ss = movingmean(Z_holding_avg,10);
t_ss = movingmean(t,10);
t_index = find(t_ss >= t_0);
[xData, yData] = prepareCurveData(t_ss(t_index), Z_holding_avg_ss(t_index) );
ft = fittype( 'a*exp(-b*x)+c', 'independent', 'x', 'dependent', 'y' );
%ft = fittype( 'exp(-a*(x+b))+c', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
%opts.Lower =  [0 0 0];
opts.StartPoint = [0.5 0.015 0.5];
[fitresult, ~] = fit( xData, yData, ft, opts );       
plot(t_ss(t_index), Z_holding_avg_ss(t_index), '.-g');
ff = plot(fitresult);
hold off
xlabel('time (s)')
ylabel('Normalized Extension (\mum)')
axis([-40 700, 0 1.2]);
fitresult
% Set the remaining axes properties
set(axes1,'FontName','Calibri','FontSize',12);
% Create axes
axes2 = axes('Parent',figure1,'Position',[0.6890625 0.11 0.2159375 0.815]);
hold(axes2,'on');

[counts,bins] = hist(Exten_end,50); %# get counts and bin locations
barh(bins,counts)
ylim([0 1.2])
set(gca,'XTick',[], 'YTick', [])

savefig(figure1,'topo1_normolized.fig');
set(gcf,'PaperPositionMode','auto')
print('topo1_normolized.png','-dpng','-r0');

f2 = figure;
Exten_end = Exten_end(~isnan(Exten_end));
 histfit(Exten_end,20)
 hold on;
 %plot the expected extension at 10 pN, such as 4106 here
 %plot([4106 4106],[0 20], '--k','LineWidth',2);
 pd = fitdist(Exten_end','Normal')
 %plot([4106 4106],[0 20], '--k','MarkerSize', 6);
 xlabel( 'Extension (\mum)');
 ylabel( 'Counts', 'Interpreter', 'none' );
 %legend('Colony#4@batch1')
 xlim([0 1.5])
 title([num2str(mean(Exten_end)) ' \pm ' num2str(std(Exten_end)) ' \mum'],'Interpreter','TEX')
savefig(f2,'extension_nomolized.fig');
set(gcf,'PaperPositionMode','auto')
print('extension_normalized.png','-dpng','-r0');

% set(gca,'FontSize',12,'FontName','Calibri');
% savefig(f2,'topo_fitted.fig');
% set(gcf,'PaperPositionMode','auto')
% print('topo_fitted.png','-dpng','-r0');

f3 = figure;
force_2 = force(index_height);
%min(height) max(height) mean(height) std(height)
h = histogram(force_2, 0:0.02:1.0)
aveForce = mean(force_2)
stdForce = std(force_2)
title({[num2str(aveForce) '\pm' num2str(stdForce) 'pN']
    ['N =' num2str(length(index_height))]},'Interpreter','TEX')
xlabel('Force (pN)')
ylabel('Count')
set(gca,'FontSize',12,'FontName','Calibri');
axis([0 1.0, 0 25])
savefig(f3,'force.fig');
set(gcf,'PaperPositionMode','auto')
print('force.png','-dpng','-r0');