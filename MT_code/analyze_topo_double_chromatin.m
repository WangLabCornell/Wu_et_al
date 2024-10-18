function analyze_topo_double_chromatin(refbeadNum)
refbeadNum = refbeadNum+1;
disp('Loading Force calibration data')
[TimeFC, ~, ~, X_FC, Y_FC, Z_FC ,~, ~] = process_MT_data('MT data_FC.txt',refbeadNum);
window = round(1/mean(diff(TimeFC)));
driftCorrX = movingmean(X_FC(:,refbeadNum),window);
driftCorrY = movingmean(Y_FC(:,refbeadNum),window);
driftCorrX = repmat(driftCorrX,1,size(Z_FC,2));
driftCorrY = repmat(driftCorrY,1,size(Z_FC,2));
X_FC = X_FC - driftCorrX;
Y_FC = Y_FC - driftCorrY;


disp('Loading Geometry calibration data')
[~, TurnGC, ~, X_GC, Y_GC, Z_GC, ~, ~] = process_MT_data('MT data_GC.txt',refbeadNum);

disp('Loading Hat curve data')
[~, TurnH, ~, ~, ~, Z_H, ~, ~] = process_MT_data('MT data_Hat.txt',refbeadNum);

disp('Loading Final hat curve data')
[~, TurnF, Box_nameF, ~, ~, Z_F, ~, ~] = process_MT_data('MT data_S.txt',refbeadNum);

disp('Loading topo1 data')
[TimeT, TurnT, ~, ~, ~, Z_T, ~, ~] = process_MT_data('MT data_topo1.txt',refbeadNum);
disp('Loading topo2 data')
[TimeT2, TurnT2, ~, ~, ~, Z_T2, ~, ~] = process_MT_data('MT data_topo2.txt',refbeadNum);

%[TimeFC, TurnFC, Box_nameFC, X_FC, Y_FC, Z_FC, L_FC, Z_piezoFC] = get_mtdata_setup2('MT data_processed_forcecalibration.txt');
%[TimeGC, TurnGC, Box_nameGC, X_GC, Y_GC, Z_GC, L_GC, Z_piezoGC] = get_mtdata_setup2('MT data_processed_geometrycorrection.txt');
%[TimeH, TurnH, Box_nameH, X_H, Y_H, Z_H, L_H, Z_piezoH] = get_mtdata_setup2('MT data_processed_Hat.txt');
%[TimeF, TurnF, Box_nameF, X_F, Y_F, Z_F, L_F, Z_piezoF] = get_mtdata_setup2('MT data_processed_finHAT.txt');

%% select tethers
% find surface
TurnF_surface = TurnF(end);
index_surface = TurnF == TurnF_surface;
surface_raw = Z_F(index_surface,:);
surface = zeros(1,size(surface_raw,2));
for i = 1:size(surface_raw,2)
    temp = surface_raw(:,i);
    temp = temp(temp > -1);
    temp = temp(temp <4);
    StDevTemp = std(temp);
    MeanTemp = mean(temp);
    temp = temp(abs(temp - MeanTemp) < 3*StDevTemp); %removes outliers that are 3 std away from the mean
    surface(i) = mean(temp);
end
%removes untracked boxes
% badTraceTracking = isnan(surface);

[TurnPartPosH, ZPartPosH,TurnPartNegH, ZPartNegH, ZStdPosH,ZStdNegH] = partHat(TurnH,Z_H);

%identify hysterysis by observing the difference between start and end of
% %hat curve.
% indexZeroTurn = TurnPartNegH == 0;
% temp = find(indexZeroTurn);
% 
% badTraceHysterysis = false(1,size(ZPartPosH,2));
% if size(temp,1) > 1
%     ZStartEnd = ZPartNegH(indexZeroTurn,:);
%     stdStartEnd = ZStdNegH(indexZeroTurn,:);
%     
%     badTraceHysterysis = (repmat(abs(diff(ZStartEnd,1,1)),2,1) - 2*stdStartEnd)>0;
%     badTraceHysterysis = badTraceHysterysis(1,:)|badTraceHysterysis(2,:);
%     
%     TurnPartNegH(temp(2)) = [];
%     ZPartNegH(temp(2),:) = [];
% end



[TurnPartPosH,I_P] = sort(TurnPartPosH);
ZPartPosH = ZPartPosH(I_P,:);
TurnNegTemp = TurnPartNegH(1:end-1);
[TurnPartNegH,I_N] = sort(TurnNegTemp);
ZPartNegTemp = ZPartNegH(1:end-1,:);
ZPartNegH = ZPartNegTemp(I_N,:);

if isequal(TurnPartPosH,TurnPartNegH)
    TurnPartH = TurnPartPosH;
else
    error('Error in hat curve partitioning : pos turn and neg turn not matching')
end

ZPartH = (ZPartPosH + ZPartNegH)./2;

% identify TC tether by hat curve 
hatHeight = max(ZPartH,[],1) - min(ZPartH,[],1);

badTraceNonTChat = hatHeight <0.1;

badTraceNonTCsurface = surface-0.1 > min(ZPartH,[],1);

badTrace_surfaceNaN = isnan(surface);

%% tethers that were TC but became non-TC during winding. 
 
 brokenTethers = find(~badTraceNonTChat & badTraceNonTCsurface);%tethers that were TC but became non-TC during winding. 
 windowBrokenTether = floor(100/(mean(diff(TimeT))));
 for i = brokenTethers
     temp = sort(movingmean(Z_T(:,i),round(windowBrokenTether/10)));
     surface(i) = mean(temp(1:windowBrokenTether));
     badTraceNonTCsurface(i) = false;
 end
 %% tethers that were TC but got detached during winding. 
 detachedTethers = ~badTraceNonTChat & badTrace_surfaceNaN;
 for i = find(detachedTethers)
     temp = movingmean(Z_T(:,i),round(windowBrokenTether/10));
      I = find(diff(Z_T(:,i))>0.3);
      if ~isempty(I)
          istart = max(1,I(end)-2-round(windowBrokenTether/10));
        surface(i) = mean(Z_T(setdiff([istart:I(end)-2],I),i));
      else
          surface(i) = min(temp);
      end
     
 end
 
surface_topo = zeros(1,size(surface_raw,2));
 for i = 1:size(surface_raw,2)
     if detachedTethers(i)
        surface_topo(i) = surface(i)  ;
     else
         [N,edges] = histcounts(Z_T(:,i),100);
         [xData, yData] = prepareCurveData(edges(1:end-1)+mean(diff(edges))/2 ,N);

        % Set up fittype and options.
        ft = fittype( 'gauss2' );
        opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
        opts.Display = 'Off';
        opts.Lower = [-Inf -Inf 0 -Inf -Inf 0];
        
        [fitresult, ~] = fit( xData, yData, ft, opts );
        
        surface_topo(i) = min(fitresult.b1,fitresult.b2);
     end
 end
 surface = min([surface ; surface_topo],[],1);
%%
badTrace_surfaceNaN = isnan(surface);


badTraceNonTC = badTraceNonTChat|badTraceNonTCsurface;


badTrace = badTraceNonTC;

badTrace = badTrace|badTrace_surfaceNaN;
if isempty(find(~badTrace,1))
    error('No good trace. Check your reference bead')
end
j = 1;
figure('Renderer', 'painters', 'Position', [623,446,560,420]);
N_goodtrace = sum(~badTrace);
for i = find(~badTrace)
    disp([num2str(j) '/' num2str(N_goodtrace)])
    scatter(TurnPartH, ZPartPosH(:,i),'MarkerFaceColor','k','MarkerEdgeColor','k')
    hold on
    scatter(TurnPartH, ZPartNegH(:,i),'MarkerFaceColor','r','MarkerEdgeColor','r')
    ylabel('Z (\mum)','FontSize',12);
    xlabel('Turn \phi','FontSize',12);
    set(gca,'FontSize',12,'FontName','Calibri');   
    GB = questdlg('Good?','Yes','No');
    if strcmp(GB,'No')
        badTrace(i) = true;
    end
    hold off
    j = j+1;
end



%remove bad traces
BoxName_sel = Box_nameF;
BoxName_sel = BoxName_sel(:,~badTrace);

Z_H_sel = Z_H;
Z_H_sel(:,badTrace) = [];
ZPartH_sel = ZPartH;
ZPartH_sel(:,badTrace) = [];

ZPartNegH_sel = ZPartNegH;
ZPartNegH_sel(:,badTrace) = [];

ZPartPosH_sel = ZPartPosH;
ZPartPosH_sel(:,badTrace) = [];


X_FC_sel = X_FC;
X_FC_sel(:,badTrace) = [];
Y_FC_sel = Y_FC;
Y_FC_sel(:,badTrace) = [];
Z_FC_sel = Z_FC;
Z_FC_sel(:,badTrace) = [];

X_GC_sel = X_GC;
X_GC_sel(:,badTrace) = [];
Y_GC_sel = Y_GC;
Y_GC_sel(:,badTrace) = [];
Z_GC_sel = Z_GC;
Z_GC_sel(:,badTrace) = [];

Z_T_sel = Z_T;
Z_T_sel(:,badTrace) = [];
Z_T2_sel = Z_T2;
Z_T2_sel(:,badTrace) = [];

Z_F_sel = Z_F;
Z_F_sel(:,badTrace) = [];

surface_sel = surface(~badTrace);

N_sel = sum(~badTrace);


%print out figures to verify selection
%f1 = figure;
%for i = 1:size(Z_H,2)
%    scatter(TurnPartH, ZPartH(:,i))
%    hold on
%   plot([min(TurnPartH) max(TurnPartH)],[surface(i) surface(i)])
%    hold off
%    if badTrace(i)
%       title('bad')
%    else
%        title('selected')
%    end
%    print([Box_nameH{i} '.png'],'-dpng','-r0');
%end
%% get geometry correction
R = 0.5;
[GCcorrection,GCcenter, ~, GCX, GCY, Rfit] = geometryCorrection(R,TurnGC,X_GC_sel,Y_GC_sel);


correction = surface_sel - GCcorrection;

ZPartH_cor = ZPartH_sel - repmat(correction,size(ZPartH_sel,1),1);
ZPartPosH_cor = ZPartPosH_sel - repmat(correction,size(ZPartH_sel,1),1);
ZPartNegH_cor = ZPartNegH_sel - repmat(correction,size(ZPartH_sel,1),1);
% ZPartPosF_cor = ZPartPosF_sel - repmat(correction,size(ZPartPosF_sel,1),1);
Z_H_cor = Z_H_sel - repmat(correction,size(Z_H_sel,1),1);
Z_FC_cor = Z_FC_sel - repmat(correction,size(Z_FC_sel,1),1);
Z_T_cor = Z_T_sel - repmat(correction,size(Z_T_sel,1),1);
Z_F_cor = Z_F_sel - repmat(correction,size(Z_F_sel,1),1);
Z_T2_cor = Z_T2_sel - repmat(correction,size(Z_T2_sel,1),1);

%% B field angle
Rfit_B = Rfit;
indexOutlier = false(1,N_sel);
for i = 1 : N_sel
   if Rfit_B(i) >0.5
       indexOutlier(i) = true;
   end
    if Rfit_B(i) <0.2
        indexOutlier(i) = true;
    end
end
[B_u,B_u_p] = getBdirection(GCX(~indexOutlier,:),GCY(~indexOutlier,:),GCcenter(~indexOutlier,:));

%% force calibration
% drift correction, projection to B field
[F_u,F_p] = forceCalibration(X_FC_sel,Y_FC_sel,Z_FC_cor,TimeFC,B_u,GCcorrection);


forces = (F_u+F_p)/2;

%% plot tether property
nucNum = zeros(1,N_sel);
isNucQualityGood = false(1,N_sel);
f3 = figure('Renderer', 'painters', 'Position', [1080 100 600 800]);
for i = 1 : N_sel
    if numel(find(isnan(ZPartH_cor(:,i))))==0
    TurnPartH_cut = TurnPartH(TurnPartH<38)
    ZPartH_cor_i = ZPartH_cor(:,i);
    ZPartH_cut = ZPartH_cor_i(TurnPartH<38)
        
    %[nucQuality,n_nuc, dn_nuc, tipSize,tipCurve] = getBraidQuality(TurnPartH,ZPartH_cor(:,i));
    [nucQuality,n_nuc, dn_nuc, tipSize,tipCurve] = getBraidQuality(TurnPartH_cut, ZPartH_cut);
    nucNum(i) = n_nuc  ;
    isNucQualityGood(i) = strcmp(nucQuality,'good');
    
    
    subplot(4,2,[1,4])
    plot(TurnPartH, ZPartH_cor(:,i),  '-ok','MarkerSize',4,'LineWidth',1.5,'MarkerFaceColor','k');
    hold on
    %plot(TurnPartNegH, ZPartNegH_cor(:,i), '-ok','MarkerSize',5,'LineWidth',2,'MarkerFaceColor','r'); 
    plot(movingmean(TurnF,100), movingmean(Z_F_cor(:,i),100),  '-o','MarkerSize',2,'LineWidth',1,'MarkerFaceColor',[0.5, 0.5, 0.5],'Color',[0.5, 0.5, 0.5]);
   
    
    ymax = max(ZPartH_cor(:,i))+0.1;
    ylim([0 max([ymax 1])]);
    titletxt = ['n_n_u_c = ' num2str(n_nuc) ' tether quality  = ' nucQuality];
    title(titletxt);
    ylabel('Z (\mum)','FontSize',12);
    xlabel('Turn','FontSize',12);
    set(gca,'FontSize',12,'FontName','Calibri');
    grid on;
    hold off
    
    
%     subplot(4,2,6)
%     scatter(X_FC(:,i), Y_FC(:,i),'MarkerFaceColor','r','MarkerEdgeColo','k')
%     hold on
%     axis equal
%     title('Force Calibration' ,'Interpreter','None','FontSize',12);
%     ylabel('Y bead position','FontSize',12);
%     xlabel('X bead position','FontSize',12);
%     set(gca,'FontSize',12,'FontName','Calibri');
%     hold off
%     
%     subplot(4,2,5)
%     plot(GCX(i,:),GCY(i,:),'s','MarkerSize',5,'LineWidth',2,'MarkerFaceColor','r');
%     hold on
%     theta = linspace(0, 2*pi,1000);
%     xcircle = GCcenter(i,1)+Rfit(i)*cos(theta);
%     ycircle = GCcenter(i,2)+Rfit(i)*sin(theta);
%     plot(xcircle,ycircle,'-','LineWidth',0.5,'Color','k');
%     title('Geometry correction' ,'Interpreter','None','FontSize',12);
%     ylabel('Y bead position','FontSize',12);
%     xlabel('X bead position','FontSize',12);
%     set(gca,'FontSize',12,'FontName','Calibri');
%     xlim([GCcenter(i,1)-0.5 GCcenter(i,1)+0.5]); ylim([GCcenter(i,2)-0.5 GCcenter(i,2)+0.5]);
    
    if abs(n_nuc-50) > 5
        figHeader = 'bad nucleosome number';
    elseif tipSize>0.07
        figHeader = 'bad tip';
    elseif tipCurve > 0.0014
        figHeader = 'bad curvature';
    else
        
        figHeader = 'good';
        
    end
        
    hold off
        
    subplot(4,2,[5 6])
    plot(TimeT,Z_T_cor(:,i));
    hold on
     xlim([0 max(TimeT)])
%      line([100 100], [-1 2],'Color','k','LineStyle','--')
%      line([470 470], [-1 2],'Color','k','LineStyle','--')
%      line([840 840], [-1 2],'Color','k','LineStyle','--')
%      line([1210 1210], [-1 2],'Color','k','LineStyle','--')
%      line([1580 1580], [-1 2],'Color','k','LineStyle','--')
     Z_T_cor_i = Z_T_cor(:,i);
     
     h0 = mean(Z_T_cor_i(TimeT <15));
     line([0 1600], [h0 h0],'Color','k','LineStyle','--')
     line([0 1600], [GCcorrection(i) GCcorrection(i)],'Color','k','LineStyle','--')
     ylim([0 min([max(Z_T_cor(:,i)) 1.5])])
    title('topo activity')
    xlabel('Turns (s)')
    ylabel('Extension (um)')
    hold off
    
    subplot(4,2,[7 8])
    plot(TimeT2,Z_T2_cor(:,i));
    hold on
     xlim([0 max(TimeT)])
      line([0 1600], [h0 h0],'Color','k','LineStyle','--')
     line([0 1600], [GCcorrection(i) GCcorrection(i)],'Color','k','LineStyle','--')
     ylim([0 min([max(Z_T_cor(:,i)) 1.5])])
    title('topo activity')
    xlabel('Turns (s)')
    ylabel('Extension (um)')
    
    hold off
    
   savefig(f3,[figHeader ' ' BoxName_sel{i} '.fig']);
   set(gcf,'PaperPositionMode','auto')
   print([figHeader ' ' BoxName_sel{i} '.png'],'-dpng','-r0');
   
    end
end
boxName = BoxName_sel;
Z_H = Z_H_cor;
Z_F = Z_F_cor;
save('Selected and Corrected traces.mat','TimeT','TurnT','Z_T_cor','TimeT2','TurnT2','Z_T2_cor','TurnH', 'Z_H','TurnF', 'Z_F','nucNum', 'isNucQualityGood', 'F_u', 'F_p','correction','badTrace','boxName','refbeadNum','GCcorrection');
%save('Selected and Corrected traces.mat','TurnH', 'Z_H', 'TurnF', 'Z_F', 'nucNum',  'isNucQualityGood', 'F_u', 'F_p','correction','badTrace','boxName','refbeadNum')
close all

end
function [F_u,F_p] = forceCalibration(X_FC,Y_FC,Z_FC,TimeFC,B_u,GCcorrection)
%% force calibration
% drift correction, projection to B field
Btheta = atan(B_u(2)/B_u(1));
B_u_p = [cos(Btheta+pi/2),sin(Btheta+pi/2)];

N_sel = size(X_FC,2);

X_FC_drift = zeros(size(X_FC));
Y_FC_drift = zeros(size(Y_FC));
F_u = zeros(1,N_sel);
F_p = zeros(1,N_sel);

kbt = 0.00408 ; %pN um
extensions = zeros(1,N_sel);
R =0.5;
heights_f = zeros(1,N_sel);
for i = 1 : N_sel
    % drift correction,
        [xData, yData] = prepareCurveData(TimeFC, X_FC(:,i));
        ft1 = fittype( 'poly1' );
        [fitresultX, ~] = fit( xData, yData, ft1 );
        
        X_FC_drift(:,i) =  X_FC(:,i) - fitresultX(TimeFC);
        
        [xData, yData] = prepareCurveData(TimeFC, Y_FC(:,i));
        ft2 = fittype( 'poly1' );
        [fitresultY, ~] = fit( xData, yData, ft2 );
        
        Y_FC_drift(:,i) =  Y_FC(:,i) - fitresultY(TimeFC);
        
        %projection to B field
        X_u_B = [X_FC_drift(:,i) Y_FC_drift(:,i)] * B_u';
        Y_u_B = [X_FC_drift(:,i) Y_FC_drift(:,i)] * B_u_p';
        
        % force calculation
        tether_height = mean(Z_FC(:,i));
        heights_f(i) = tether_height;
        varX = var(X_u_B);
        varY = var(Y_u_B);
        F_u(i) = kbt *tether_height/varX;
        F_p(i) = kbt *(tether_height+R-GCcorrection(i))/varY;
        extensions(i) = tether_height*1000;
        
end
end
function [Times, Turn, Box_name, X_bead, Y_bead, Z_bead_d_b2b ,L, Z_piezo, isCorrected] = process_MT_data(filename,refbeadNum,zscanName)

[Times, Turn, Box_name, X_bead, Y_bead, Z_bead,L, Z_piezo] = get_mtdata_raw(filename);
%pick a bead
%refbeadNum = 2;
%Z_bead_fs = Z_bead;
Z_ref_bead = Z_bead(:,refbeadNum);
dt = mean(diff(Times));
windowSize = round(1/dt);

Z_bead_b2b = Z_bead;
isCorrected = false(1,size(Z_bead,2));
if nargin == 3
    parfor i = 1:size(Z_bead,2)
        [isCorrected(i),Z_bead_b2b(:,i)] = correctB2Bvar(Z_bead(:,i), i, refbeadNum, zscanName);
    end
else
    Z_bead_b2b = Z_bead;
end

Z_bead_d_b2b = correctZdrift(Z_bead_b2b, Z_ref_bead, windowSize);


Z_bead_d_b2b = Z_bead_d_b2b*1.2*1.33;% For MT1

Z_bead_d_b2b = removeTrackingError(Z_bead_d_b2b,-1);



%% drift correct

%% selct beads that  1. do not have tracking error 2. wound to surface
%% use raw data to correct b2b variation
%% drift correct again to get final result. 




end
function [Times, Turn, Box_name, X_bead, Y_bead, Z_bead,L, Z_piezo] = get_mtdata_raw(filename)
% use raw data without correction

A = importdata(filename);


    Times = A.data(:,1)/1000;           %s 
    Turn = A.data(:, 3)/32; 
    isData = false;
    colNum = 0;
    while ~isData
        colNum = colNum+1;
        temp = A.colheaders(colNum);
        temp = temp{1};
        isData = strcmp(temp(1:3),'Box');
    end
        
    Box_data = A.data(:,colNum:end);
    Box_name = A.colheaders(colNum:3:end);
    L = length(Box_name);
    X_bead = Box_data(:,1:3:end) * 0.1465 *2;
    Y_bead = Box_data(:,2:3:end) * 0.1465 *2;
    Z_bead = Box_data(:,3:3:end);
    Z_piezo = A.data(:,4);    

end
function zBeadCorr = correctZdrift(zBead, zrefBead, windowSize) 

Zcorrection = movingmean(zrefBead,windowSize);
zBeadCorr = zBead - repmat(Zcorrection,1,size(zBead,2));
end
function result = removeTrackingError(Zbead,threshold)
    indexTE = Zbead < threshold | isnan(Zbead);
    
    patches = indexTE(2:end,:) - indexTE(1:end-1,:);
    patches = [zeros(1,size(Zbead,2)); patches];
    result = Zbead;
    
    for k = 1:size(Zbead,2)
        
    starts = find(patches(:,k) == 1);
    ends = find(patches(:,k) == -1)-1;
    
    numError = numel(starts);
    
    if ~isequal(size(starts),size(ends)) % when the experiment ends with a tracking error 
        ends = [ends; size(Zbead,1)];
        numError = numError-1;
        if numel(starts)>0
        if starts(end) == ends(end)
             j = starts(end);
            result(j,k) = result(j-1,k);
        else
            v = repmat(Zbead(starts(end)-1,k),ends(end)-starts(end)+1,1);
            result(starts(end):ends(end),k) = v;
        end
        end
    end
    
        
    
    if numError > 0
        for i = 1 :numError
            if starts(i) == ends(i)
                j = starts(i);
                result(j,k) = (result(j-1,k) + result(j+1,k))*0.5;
            else
                v = linspace(Zbead(starts(i)-1,k), Zbead(ends(i)+1,k),ends(i)-starts(i)+1);
                result(starts(i):ends(i),k) = v;
            end
        end
    end
        
    end

end
function [TurnPositive, ZPositive, TurnNegative, ZNegative,StdPositive,StdNegative] = partHat(Turn,ZBead)
%% partition hat curve data to discrete turns
% assumes that hat curve data was taken in following order : 0 --> Min
% negative turns --> max Positive turns -->0
k=1;
%time_pat{1} = Time(1);
turn_pat{1} = Turn(1);
index_pat{1} = 1;
for i = 1 : length(Turn)-1
        if abs(Turn(i+1)-turn_pat{k}(end))<0.01
            %time_pat{k} = [time_pat{k} Time(i+1)];   
            turn_pat{k} = [turn_pat{k} Turn(i+1)];   
            index_pat{k} = [index_pat{k} i+1];   
        elseif abs(Turn(i+1) - turn_pat{k}(end)) == 1
            k = k+1;
            %time_pat{k} = Time(i+1);
            turn_pat{k} = Turn(i+1);
            index_pat{k} = i+1;
        end
end    
TurnMean  = zeros(length(turn_pat),1);
ZMean  = zeros(length(turn_pat),size(ZBead,2));
ZStd = ZMean;
for i = 1 : size(TurnMean,1)
        TurnMean(i) = mean(turn_pat{i}); 
        ZMean(i,:) = mean(ZBead(index_pat{i},:),1);
        ZStd(i,:) = std(ZBead(index_pat{i},:),1);
end

index_negativeSteps = diff(TurnMean)<0;
index_positiveSteps = ~index_negativeSteps;

temp = find(index_positiveSteps);
index_negativeSteps(temp(1)) = true;
index_negativeSteps(end+1) = true;
index_positiveSteps(temp(end)+1) = true;


TurnPositive = TurnMean(index_positiveSteps);
ZPositive = ZMean(index_positiveSteps,:);
StdPositive = ZStd(index_positiveSteps,:);
TurnNegative = TurnMean(index_negativeSteps); 
ZNegative = ZMean(index_negativeSteps,:);
StdNegative = ZStd(index_negativeSteps,:);

end
function [correction,center, GCangles, GCX, GCY,Rfit] = geometryCorrection(R,TurnGC,XBead,YBead)
    %ver 2020/12/23
     Ntraces = size(XBead,2);
    thetas = 0:0.125:0.875;
    
    GCangles = thetas*2*pi;
    Nangles = numel(thetas);
    GCX = zeros(Nangles, Ntraces);
    GCY = zeros(Nangles, Ntraces);
    
    %% correct drift
    
    finAngle = find(TurnGC == 0.875);
    XBead = XBead(1:finAngle(end),:);
    YBead = YBead(1:finAngle(end),:);
    theta = TurnGC(1:finAngle(end));   
    XBead = XBead(ismember(theta,thetas),:);
    YBead = YBead(ismember(theta,thetas),:);
    theta = theta(ismember(theta,thetas));
    
    XdriftAll = zeros(size(theta));
    XdriftPrev = 0;
    YdriftAll = zeros(size(theta));
    YdriftPrev = 0;
    
    for i = 1:Nangles
        i
        thisTheta = theta == thetas(i);
        %Xdrift
        temp = movingmean(mean(XBead(thisTheta,:),2),min(numel(mean(XBead(thisTheta,:),2)),50));
        XdriftThis = temp - temp(1)+XdriftPrev(end);
        XdriftAll(thisTheta) = XdriftThis; 
        XdriftPrev = XdriftThis;
        %Ydrift
        temp = movingmean(mean(YBead(thisTheta,:),2),50);
        YdriftThis = temp - temp(1)+YdriftPrev(end);
        YdriftAll(thisTheta) = YdriftThis; 
        YdriftPrev = YdriftThis;
    end
    
    XBead_cor = XBead - repmat(XdriftAll,1,Ntraces);
    YBead_cor = YBead - repmat(YdriftAll,1,Ntraces);
    %%
    
    for i = 1:Nangles
    
        thisTheta = theta == thetas(i);
        GCX(i,:) = mean(XBead_cor(thisTheta,:),1);
        GCY(i,:) = mean(YBead_cor(thisTheta,:),1);
    end
    
    center = zeros(Ntraces,2);
    Rfit = zeros(1,Ntraces);
    
    for i = 1 : Ntraces
       a=[GCX(:,i) GCY(:,i) ones(size(GCX(:,i)))]\(-(GCX(:,i).^2+GCY(:,i).^2));
       center(i,:) = [-.5*a(1), -.5*a(2)];
       Rfit(i)  =  sqrt((a(1)^2+a(2)^2)/4-a(3));
    end
    correction = R -  heaviside(R-Rfit) .* sqrt(R^2 - Rfit.^2);
    GCX = GCX';
    GCY = GCY';

end
function [B_u,B_u_p] = getBdirection(GCX,GCY,GCcenter)
%% B field angle
BeadDir = [GCX(:,1) GCY(:,1)] - GCcenter;

BeadAngle = atan((BeadDir(:,2)+0.1)./(BeadDir(:,1)));
index = BeadAngle<0;
BeadAngle(index) = BeadAngle(index);
Btheta = mean(BeadAngle);
slope = tan(Btheta);

B_u = [cos(Btheta),sin(Btheta)]; %unit vector in the direction of B field
B_u_p = [cos(Btheta+pi/2),sin(Btheta+pi/2)]; %unit vector in the perpendicular direction of B field

figure;
hold all
scatter( BeadDir(:,1),   BeadDir(:,2))
scatter(0,0)
%for i = 1 : size(BeadAngle,1)
%    scatter( BeadAngle(i,1),   BeadAngle(i,2))
%    plot([0 BeadAngle(i,1)], [0, BeadAngle(i,2)], 'Color' ,'r','Marker','v');
%end
xrange = -0.5:0.01:0.5;
yfit = xrange*slope;
plot(xrange,yfit)
xlim([-0.5 0.5])
ylim([-0.5 0.5])
hold off
print(['Bead Angle.png'],'-dpng','-r0');
end

function [nucQuality,n_nuc, dn_nuc, tipSize,curvature] = getBraidQuality(Turn_ev,Z_ev)
ft = fittype( 'poly2' );
indexMax = find(Z_ev == max(Z_ev));
indexMax = indexMax(1);
Turn_max = Turn_ev(indexMax);
maxHeight = Z_ev(indexMax);

index_fitRange = find(Turn_ev > Turn_max-10 & Turn_ev < Turn_max+10);
Turn_fit = Turn_ev(index_fitRange);
Z_fit = Z_ev(index_fitRange);
[xData, yData] = prepareCurveData(Turn_fit, Z_fit);
[fitresult, ~] = fit( xData, yData, ft );

tipSize = (maxHeight - fitresult(Turn_max))/maxHeight;
curvature=abs(fitresult.p1);
dL = 200;
%TUngs 'universal equation'
%No = 12644; 

%load('x_f_WLC.mat')
%n_nuc =  (maxHeight*1000 - No * 0.338 * interp1(Force_WLC,L_r, force/2))/(7.1*force/2+1.23 - 197 * 0.338 * interp1(Force_WLC,L_r, force/2));
%dn_nuc = dL/abs(7.1*force/2+1.23 - 197 * 0.338 * interp1(Force_WLC,L_r, force/2) );        


slope = -0.0371;
intercept = 2.8031;
n_nuc =  (maxHeight - intercept)/slope;
dn_nuc = abs(dL/slope); 

isNucNumGood = abs(n_nuc-50) <= 5;
isTipGood = tipSize <0.07;
isCurveGood = abs(fitresult.p1) < 0.0014;

if isNucNumGood && isTipGood && isCurveGood
    nucQuality = 'good';
else
    nucQuality = 'bad';
end
end
