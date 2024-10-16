
% load data
clear
load('Selected and Corrected traces.mat')
%trace{1}.sign_winding

% Select good Nuc traces within 50+/-5 and Extract parameters 
index = 1;
Z_T_cor_good = [];
Z_T2_cor_good = [];
force = [];
height = [];
for i = 1:numel(isNucQualityGood)
    if isNucQualityGood(i)
        Z_T_cor_i = Z_T_cor(:,i);
        temp_height_0 = mean(Z_T_cor_i(TimeT<=10));
        Z_T_cor_good(:,index) = Z_T_cor(:,i)/temp_height_0;%normolized height
        Z_T2_cor_good(:,index) = Z_T2_cor(:,i)/temp_height_0;%normolized height
        force(index) = F_u(i);
        %height(index) = extensions(i);
        index = index + 1;
    end
end

Exten_end = zeros(numel(force),0);
% Create figure
figure1 = figure;
% Create axes
axes1 = axes('Parent',figure1,...
    'Position',[0.13 0.11 0.542916666666667 0.815]);
hold(axes1,'on');

t1 = TimeT - 40;
for i = 1:numel(force)
    
      dt = mean(diff(TimeT));
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
    Z_T_cor_good_temp = sgolayfilt(Z_T_cor_good(:,i),order,window2);
    Exten_end(i) = mean(Z_T_cor_good_temp(end-0.1*100:end));
    plot1 = plot(t1, Z_T_cor_good_temp,'Color',[0.800000011920929 0.800000011920929 0.800000011920929]);
    %plot(time_holding_1(:,1),z_holding_1(:,i))
end
Z_T_cor_good_avg = nanmean(Z_T_cor_good,2);
plot(t1, Z_T_cor_good_avg)
%line([t(turn_index(1)) t(turn_index(1))], [1.5 4],'LineStyle','--','LineWidth',1, 'Color','k');
line([ -40 700], [1 1],'LineStyle','--','LineWidth',1, 'Color','k');
%turn_index = find(TurnT == 70);
%t_0 = TimeT(turn_index(1))+10;

% f2 = figure;
% plot(t1, Z_T_cor_good_avg)
% %turn_index = find(TurnT == 70);
 t_0 = 45;
% t_index = t(turn_index(1)).
%hold on;
Z_T_cor_good_avg_ss = movingmean(Z_T_cor_good_avg,10);
t_ss = movingmean(t1,10);
t_index = find(t_ss >= t_0);
[xData, yData] = prepareCurveData(t_ss(t_index), Z_T_cor_good_avg_ss(t_index) );
ft = fittype( 'a*exp(-b*x)+c', 'independent', 'x', 'dependent', 'y' );
%ft = fittype( 'exp(-a*(x+b))+c', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
%opts.Lower =  [0 0 0];
opts.StartPoint = [1 0.005 2];
[fitresult, ~] = fit( xData, yData, ft, opts );       
plot(t_ss(t_index), Z_T_cor_good_avg_ss(t_index), '.-g');
ff = plot(fitresult);
%line([ -40 700], [mean(height) mean(height)],'LineStyle','--','LineWidth',1, 'Color','k');
%line([t1(turn_index(1)) t1(turn_index(1))], [0 1.8],'LineStyle','--','LineWidth',1, 'Color','k');
hold off
xlabel('time (s)')
ylabel('Extension (\mum)')
axis([-40 700, 0 1.2]);
fitresult
% Set the remaining axes properties
set(axes1,'FontName','Calibri','FontSize',12);
% Create axes
axes2 = axes('Parent',figure1,'Position',[0.6890625 0.11 0.2159375 0.815]);
hold(axes2,'on');

[counts,bins] = hist(Exten_end,20); %# get counts and bin locations
barh(bins,counts)
ylim([0 1.2])
set(gca,'XTick',[], 'YTick', [])

set(gca,'FontSize',12,'FontName','Calibri');
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



f3 = figure;
%min(height) max(height) mean(height) std(height)
h = histogram(force, 0.3:0.02:0.7)
aveForce = mean(force)
stdForce = std(force)
title({[num2str(aveForce) '\pm' num2str(stdForce) 'pN']
    ['N =' num2str(length(force))]},'Interpreter','TEX')
xlabel('Force (pN)')
ylabel('Count')
set(gca,'FontSize',12,'FontName','Calibri');
xlim([0.3 1]);
savefig(f3,'force.fig');
set(gcf,'PaperPositionMode','auto')
print('force.png','-dpng','-r0');