function analyze_topo_dna_condensation
close all

%surface calibration before topo 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 1: Read data and input parameter
R_bead = 0.5;
kbt = 4.08; % pN.nm
rise = 0.338; % nm
Lo13p7 = 13700;    
f_e = 0.5; %pN

load('x_f_WLC.mat');
[Times1, Turn1, Box_name1, X_bead1, Y_bead1, Z_bead1, L1, Z_piezo1] = get_mtdata4('MT data with MT DAQ data_forcecalibration.txt');
[Timesg, Turng, Box_nameg, X_beadg, Y_beadg, Z_beadg, Lg, Z_piezog] = get_mtdata4('MT data with MT DAQ data_geometrycorrection.txt');
[Times2, Turn2, Box_name2, X_bead2, Y_bead2, Z_bead2, L2, Z_piezo2] = get_mtdata4('MT data with MT DAQ data_initialhat.txt');

[Times4, Turn4, Box_name4, X_bead4, Y_bead4, Z_bead4, L4, Z_piezo4] = get_mtdata4('MT data with MT DAQ data_finalhat.txt');
%[Timesf, Turnf, Box_namef, X_beadf, Y_beadf, Z_beadf, Lf, Z_piezof] = get_mtdata4('MT data with MT DAQ data_surf.txt');
% surface_index = find(Turn4 == max(Turn4));
% Times4 = Times4(1:surface_index(end));
% Turn4 = Turn4(1:surface_index(end));
% Z_bead4 = Z_bead4(1:surface_index(end),:);
% X_bead4 = X_bead4(1:surface_index(end),:);
% Y_bead4 = Y_bead4(1:surface_index(end),:);
% Z_piezo4 = Z_piezo4(1:surface_index(end),:);
%[Times5, Turn5, Box_name5, X_bead5, Y_bead5, Z_bead5, L5, Z_piezo5] = get_mtdata3('MT data_processed_finalHAT 2.txt');

%Turn4 = [Turn4; Turn5];
%X_bead4 = [X_bead4; X_bead5];
%Y_bead4 = [Y_bead4; Y_bead5];
%Z_bead4 = [Z_bead4; Z_bead5];
%Z_piezo4 = [Z_piezo4; Z_piezo5];
%Times4 = [Times4; Times5 + Times4(end)];

% first winding
[Times3, Turn3, Box_name3, X_bead3, Y_bead3, Z_bead3, L3, Z_piezo3] = get_mtdata4('MT data with MT DAQ data_topo1.txt');

% second winding
[Times5, Turn5, Box_name5, X_bead5, Y_bead5, Z_bead5, L5, Z_piezo5] = get_mtdata4('MT data with MT DAQ data_topo2.txt');
% third winding
%[Times6, Turn6, Box_name6, X_bead6, Y_bead6, Z_bead6, L6, Z_piezo6] = get_mtdata3('MT data_processed_topo3.txt');

% check local hat
%[Times7, Turn7, Box_name7, X_bead7, Y_bead7, Z_bead7, L7, Z_piezo7] = get_mtdata3('MT data_processed_checkafterhat.txt');

j = 0;
for i = 1 : L1
    if strcmp(Box_name1{i}, 'Box222.1 X (px)') %|| strcmp(Box_name1{i}, 'Box249.1 X (px)') || strcmp(Box_name1{i}, 'Box251.1 X (px)')  %7
        j = j+1;
        index_stuck(j) = i;
    end
end
index_stuck

Z_beadg = Z_beadg - repmat(mean(Z_beadg(:, index_stuck),2),1,Lg);
Z_bead1 = Z_bead1 - repmat(mean(Z_bead1(:, index_stuck),2),1,L1);
Z_bead2 = Z_bead2 - repmat(mean(Z_bead2(:, index_stuck),2),1,L2);
Z_bead3 = Z_bead3 - repmat(mean(Z_bead3(:, index_stuck),2),1,L3);
Z_bead4 = Z_bead4 - repmat(mean(Z_bead4(:, index_stuck),2),1,L4);
Z_bead5 = Z_bead5 - repmat(mean(Z_bead5(:, index_stuck),2),1,L5);
%Z_beadf = Z_beadf - repmat(mean(Z_beadf(:, index_stuck),2),1,Lf);
%% Part 2: Select TC tethers, surface finding, geometry correction, find B-field

%% calculate z_variance during second twisting to identify non-tc tether
z_variance = mean(Z_bead2.^2,1) - mean(Z_bead2,1).^2;       
%% find tether with smaller variance and assign them as tether with chromatin
index_chromatin = mean(Z_bead2,1)>-2 & mean(Z_bead3,1)>-2 & mean(Z_bead4,1)>-2 & z_variance>0.04  ;%& z_variance<0.0850*(0.8*12644*0.338/1000)^2 ;
L = length(find(index_chromatin));
Box_namec1 = Box_name1(index_chromatin); Box_namec2 = Box_name2(index_chromatin); Box_namec3 = Box_name3(index_chromatin); Box_namec4 = Box_name4(index_chromatin);
Box_namec5 = Box_name5(index_chromatin); %Box_namec6 = Box_name6(index_chromatin); Box_namec7 = Box_name7(index_chromatin); 

X_beadc1 = X_bead1(:,index_chromatin); X_beadc2 = X_bead2(:,index_chromatin); X_beadc3 = X_bead3(:,index_chromatin); X_beadc4 = X_bead4(:,index_chromatin);
Y_beadc1 = Y_bead1(:,index_chromatin); Y_beadc2 = Y_bead2(:,index_chromatin); Y_beadc3 = Y_bead3(:,index_chromatin); Y_beadc4 = Y_bead4(:,index_chromatin);
Z_beadc1 = Z_bead1(:,index_chromatin); Z_beadc2 = Z_bead2(:,index_chromatin); Z_beadc3 = Z_bead3(:,index_chromatin); Z_beadc4 = Z_bead4(:,index_chromatin);

X_beadc5 = X_bead5(:,index_chromatin); %X_beadcf = X_beadf(:,index_chromatin); %X_beadc7 = X_bead7(:,index_chromatin); 
Y_beadc5 = Y_bead5(:,index_chromatin); %Y_beadcf = Y_beadf(:,index_chromatin); %Y_beadc7 = Y_bead7(:,index_chromatin); 
Z_beadc5 = Z_bead5(:,index_chromatin); %Z_beadcf = Z_beadf(:,index_chromatin); %Z_beadc7 = Z_bead7(:,index_chromatin); 


%% find the surface, applied geometry correction from the stuck chromatin
surface = mean(Z_beadc4(Times4>Times4(end)-5-5 & Times4<Times4(end)-5,:),1); % choose turn where the chromatin were clearly stuck
 %surface_index = Turnf == 150;
 %surface = mean(Z_beadcf(surface_index,:),1); % choose turn where the chromatin were clearly stuck
Z_beadc1 = Z_beadc1 - repmat(surface,size(Z_beadc1,1),1);
Z_beadc2 = Z_beadc2 - repmat(surface,size(Z_beadc2,1),1);
Z_beadc3 = Z_beadc3 - repmat(surface,size(Z_beadc3,1),1);
Z_beadc4 = Z_beadc4 - repmat(surface,size(Z_beadc4,1),1);
Z_beadc5 = Z_beadc5 - repmat(surface,size(Z_beadc5,1),1);
%Z_beadc6 = Z_beadc6 - repmat(surface,size(Z_beadc6,1),1);
%Z_beadc7 = Z_beadc7 - repmat(surface,size(Z_beadc7,1),1);


Z_beadc2 = z_correct(Times3, Z_beadc2, L);
Z_beadc3 = z_correct(Times3, Z_beadc3, L);
Z_beadc4 = z_correct(Times3, Z_beadc4, L);
Z_beadc5 = z_correct(Times5, Z_beadc5, L);

%% find anchoring pos
pos_tether = zeros(L,2); 
%% dealing with drift
for j = 1 : L
        [xData, yData] = prepareCurveData( Times1, X_beadc1(:,j));
        ft1 = fittype( 'poly1' );
        [fitresult1, ~] = fit( xData, yData, ft1 );
        [xData, yData] = prepareCurveData( Times1, Y_beadc1(:,j));
        ft2 = fittype( 'poly1' );
        [fitresult2, ~] = fit( xData, yData, ft2 );
        pos_tether(j, :) =  [fitresult1(Times1(end)), fitresult2(Times1(end))];
end
               
%% find relative anchoring point
% partition the data according to turn number for the twisting in step 2 (Geometry correction)
k=1;
time_patg={}; turn_patg={};index_patg={};
time_patg{1} = Timesg(1);
turn_patg{1} = Turng(1);
index_patg{1} = 1;
for i = 1 : length(Turng)-1
        if abs(Turng(i+1)- turn_patg{k}(end))<0.01
            time_patg{k} = [time_patg{k} Timesg(i+1)];   
            turn_patg{k} = [turn_patg{k} Turng(i+1)];   
            index_patg{k} = [index_patg{k} i+1];   
        else
            if Turng(i+1)-turn_patg{k}(end) >= 0.12
                k = k+1;
                time_patg{k} = Timesg(i+1);
                turn_patg{k} = Turng(i+1);
                index_patg{k} = i+1;
            end
        end
end
    
K = length(index_patg);
colorK = rand(K,3);    
n1 = length(time_patg);
Turng_mean  = zeros(1,n1);
X_beadcg_mean  = [];
Y_beadcg_mean  = []; 
X_beadcg = X_beadg(:,index_chromatin);
Y_beadcg = Y_beadg(:,index_chromatin);
for i = 1 : n1
        Turng_mean(i) = mean(turn_patg{i}); 
        X_beadcg_mean(i,:) = mean(X_beadcg(index_patg{i},:),1);
        Y_beadcg_mean(i,:) = mean(Y_beadcg(index_patg{i},:),1);
end

%% fit circle to find center_rotation usign Kasa algorithm
% source: https://www.mathworks.com/matlabcentral/fileexchange/5557-circle-fit
% more accurate algorithm: Pratt and Taubin circle fit
center_rotation = zeros(size(X_beadcg_mean,2), 2);
Rfit = zeros(size(X_beadcg_mean,2),1);
for i = 1 : size(X_beadcg_mean,2)
       a=[X_beadcg_mean(:,i) Y_beadcg_mean(:,i) ones(size(X_beadcg_mean(:,i)))]\(-(X_beadcg_mean(:,i).^2+Y_beadcg_mean(:,i).^2));
       center_rotation(i,:) = [-.5*a(1), -.5*a(2)];
       Rfit(i)  =  sqrt((a(1)^2+a(2)^2)/4-a(3));
end
pos_anchoring_relative_center = center_rotation - pos_tether;        % vector pointing from bead center to the anchoring point
radius = sqrt(pos_anchoring_relative_center(:,1).^2 + pos_anchoring_relative_center(:,2).^2);
height_correction = R_bead -  heaviside(R_bead-radius) .* sqrt(R_bead^2 - radius.^2);
%% correction to height 
Z_beadc1 = Z_beadc1+ repmat(height_correction',size(Z_beadc1,1),1);
Z_beadc2 = Z_beadc2+ repmat(height_correction',size(Z_beadc2,1),1);
Z_beadc3 = Z_beadc3+ repmat(height_correction',size(Z_beadc3,1),1);
Z_beadc4 = Z_beadc4+ repmat(height_correction',size(Z_beadc4,1),1);
Z_beadc5 = Z_beadc5+ repmat(height_correction',size(Z_beadc5,1),1);
%Z_beadc6 = Z_beadc6+ repmat(height_correction',size(Z_beadc6,1),1);
%Z_beadc7 = Z_beadc7+ repmat(height_correction',size(Z_beadc7,1),1);
%% Obtain the real-time hat
%index_qw = find( Turn3 == -40);
%t_r = Times3(index_qw(end));  % time when the quick winding starts
%Turn_r = Turn3(Times3 >= t_r); X_bead_r = X_bead3(Times3 >= t_r,:); Y_bead_r = Y_bead3(Times3 >= t_r,:); Z_bead_r = Z_bead3(Times3 >= t_r,:); L_r = L3;  Z_piezo_r = Z_piezo3(Times3 >= t_r,:); Box_name_r = Box_name3;
%X_beadc_r = X_beadc3(Times3 >= t_r,:); Y_beadc_r = Y_beadc3(Times3 >= t_r,:); Z_beadc_r = Z_beadc3(Times3 >= t_r,:); Times_r = Times3(Times3 >= t_r);

%% typical topo reaction region
%Turn3 = Turn3(Times3 < t_r); X_bead3 = X_bead3(Times3 < t_r,:); Y_bead3 = Y_bead3(Times3 < t_r,:); Z_bead3 = Z_bead3(Times3 < t_r,:); Z_piezo3 = Z_piezo3(Times3 < t_r,:);
%X_beadc3 = X_beadc3(Times3 < t_r,:); Y_beadc3 = Y_beadc3(Times3 < t_r,:); Z_beadc3 = Z_beadc3(Times3 < t_r,:); Times3 = Times3(Times3 < t_r);

%% finding B field
x_anchor = pos_anchoring_relative_center(:,1);
y_anchor = pos_anchoring_relative_center(:,2);
[xData, yData] = prepareCurveData( x_anchor, y_anchor );
ft = fittype( 'poly1' );
opts = fitoptions( 'Method', 'LinearLeastSquares' );
opts.Robust = 'Bisquare';
[fitresult_Bfield, gof] = fit( xData, yData, ft, opts );
uB = [1, fitresult_Bfield.p1]/norm([1, fitresult_Bfield.p1]); % direction of B field
  
f1 = figure;
plot(0, 0,'xk','MarkerSize',10, 'LineWidth',5);
hold all
xas = -0.8:0.01:0.8;
yas = fitresult_Bfield(xas);
plot(xas, yas, '-k','LineWidth',1.5);
hold all
for i = 1 : L
        plot([0 pos_anchoring_relative_center(i,1)], [0, pos_anchoring_relative_center(i,2)],'LineWidth',1.5, 'Color', 'r', 'Marker','v', 'MarkerFaceColor','w');
end
xlim([-0.5 0.5]);
ylim([-0.5 0.5]);
axis equal
xlabel('\DeltaX(\mum)');
ylabel('\DeltaY(\mum)');
set(gca,'FontSize',12,'FontName','Calibri');


    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 4: Average before and after data during the pulse
%% average data during the pulse
%% partition the before topo data

% before hat
[Turn2_mean_1way, Z_beadc2_mean_1way, Turn2_mean_1way_negative, Z_beadc2_mean_1way_negative] = average_at_pause(Times2, Turn2, Z_beadc2);
% final hat
[Turn4_mean_1way, Z_beadc4_mean_1way, Turn4_mean_1way_negative, Z_beadc4_mean_1way_negative] = average_at_pause(Times4, Turn4, Z_beadc4);

% check after hat
%[Turn7_mean_1way, Z_beadc7_mean_1way, Turn7_mean_1way_negative, Z_beadc7_mean_1way_negative] = average_at_pause(Times7, Turn7, Z_beadc7);
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 5: Run through each trace, manually select good trace, and identify HAT shape and save statistics
%% run through each trace and save data
uB_var1 = zeros(1,L);
uB_var2 = zeros(1,L);
tv = 60;
%change here
% Turn3_600 = Turn3(Times3 < 610);
 d_wound1 = Turn3(end);
 %Turn3_index = find(Turn3 == d_wound1);
% Z_beadc3_wait = Z_beadc3(1:Turn3_index(end), :);

d_wound2 = Turn5(end);
%d_wound3 = Turn6(end);

trace = cell(1,L);
for j = 1:L
   %if j == 13 || j == 41 || j == 51 || j == 53 || j == 63 || j == 108 ||...
   % j == 115 || j == 125 || j == 137 || j == 138 || j == 153
   %   j = j+1;
   % else 
   %   j = j;
   %end
    if 1%strcmp(Box_namec1{j}, 'Box27.0 X (px)') 
        trace{j}.name = Box_namec1{j};
        ftemp = figure;
        hold all   
        
        line([0 0],[ -2 15],'LineStyle','--','LineWidth',2, 'Color','k');
        line([d_wound1 d_wound1],[ -2 15],'LineStyle','--','LineWidth',2, 'Color',[0.87, 0.49, 0]);
        %line([d_wound1 d_wound1],[ -2 15],'LineStyle','--','LineWidth',2, 'Color','b');
        %line([d_wound1 d_wound1],[ -2 15],'LineStyle','--','LineWidth',2, 'Color',[0, 0.5, 0]);
                
        %plot(Turn4_mean_1way_negative, Z_beadc4_mean_1way_negative(:,j), '-or','MarkerSize',5,'LineWidth',1.5,'Color',[0.5, 0.5, 0.5]);
        
        plot(movingmean(Turn3,10), movingmean(Z_beadc3(:,j),10), '-','LineWidth',2,'Color',[0.87, 0.49, 0]);    
        %plot(movingmean(Turn5,10), movingmean(Z_beadc5(:,j),10), '-','LineWidth',2,'Color','b');    
        %plot(movingmean(Turn6,10), movingmean(Z_beadc6(:,j),10), '-','LineWidth',2,'Color',[0, 0.5, 0]);    
        
        plot(Turn2_mean_1way, Z_beadc2_mean_1way(:,j), '-ok','MarkerSize',5,'LineWidth',2,'MarkerFaceColor','k');
        
        plot(Turn4_mean_1way, Z_beadc4_mean_1way(:,j), '-or','MarkerSize',5,'LineWidth',1.5,'Color',[1, 0.6, 0.78]);
       
        %plot(Turn3(Turn3_index(end):end), Z_beadc3(Turn3_index(end):end,j), '.','Color',[0.87, 0.49, 0]);    
        
        title(Box_namec1{j} ,'Interpreter','None','FontSize',12);
        ylabel('Z (\mum)','FontSize',12);
        xlabel('Turn','FontSize',12);
        set(gca,'FontSize',12,'FontName','Calibri');
        hmm = 2*max(Z_beadc2_mean_1way(:,j));
        xlim([-60 200]); ylim([min([-0.5 hmm+0.2]), max([-0.5 hmm+0.2]) ]);
        grid on;

        qstring = 'Does the trace look good?';
        choice = questdlg(qstring,'Question?','Yes','No','Yes');

        if strcmp(choice, 'Yes') 
            close(ftemp);

            %% Force calculation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
            [height_c_before, F_from_var_length, F_from_var_length_p] = force_calculation(Times1, X_beadc1(:, j), Y_beadc1(:,j), Z_beadc1(:,j), uB, height_correction(j), tv);
            %trace{j}.height_c_before = height_c_before;

            % initial HAT
            trace{j}.turn_m_before_positive = Turn2_mean_1way;
            trace{j}.Z_m_before_positive = Z_beadc2_mean_1way(:,j);
            trace{j}.turn_m_before_negative = Turn2_mean_1way_negative;
            trace{j}.Z_m_before_negative = Z_beadc2_mean_1way_negative(:,j);
            trace{j}.height_c_before = max(Z_beadc2_mean_1way(:,j)) ;

            % after hat
            trace{j}.turn_m_after_positive = Turn4_mean_1way;
            trace{j}.Z_m_after_positive = Z_beadc4_mean_1way(:,j);
            trace{j}.turn_m_after_negative = Turn4_mean_1way_negative;
            trace{j}.Z_m_after_negative = Z_beadc4_mean_1way_negative(:,j);
            trace{j}.height_c_after = max(Z_beadc4_mean_1way(:,j)) ;

            % force 
            trace{j}.radius = radius(j);
            trace{j}.F_from_var_length = F_from_var_length;
            trace{j}.F_from_var_length_p = F_from_var_length_p;

             % topo activity, first winding
            trace{j}.time_holding1 = Times3; time_holding1 = Times3;
             trace{j}.turn_holding1 = Turn3; turn_holding1 = Turn3;
            trace{j}.z_holding1 = Z_beadc3(:,j); z_holding1 = Z_beadc3(:,j);

             % topo activity, second winding
             trace{j}.time_holding2 = Times5; time_holding2 = Times5;
            trace{j}.turn_holding2 = Turn5; turn_holding2 = Turn5;
             trace{j}.z_holding2 = Z_beadc5(:,j); z_holding2 = Z_beadc5(:,j);
           

            % topo activity, third winding
            %trace{j}.time_holding3 = Times6; time_holding3 = Times6;
            %trace{j}.turn_holding3 = Turn6; turn_holding3 = Turn6;
            %trace{j}.z_holding3 = Z_beadc6(:,j); z_holding3 = Z_beadc6(:,j);

            % check local after hat
            %trace{j}.turn_m_checkafter_positive = Turn7_mean_1way;
            %trace{j}.Z_m_checkafter_positive = Z_beadc7_mean_1way(:,j);
            %trace{j}.turn_m_checkafter_negative = Turn7_mean_1way_negative;
            %trace{j}.Z_m_checkafter_negative = Z_beadc7_mean_1way_negative(:,j);
            %trace{j}.height_c_checkafter = 1/2 * (max(Z_beadc7_mean_1way(:,j)) + max(Z_beadc7_mean_1way_negative(:,j)));
            

            %% Fit the hat, deduce height-turn state conversion on the right and the left sides
            [paraf, height_0, turn_0, ttempt, ztempt,fitresult_height_turn_right, fitresult_height_turn_left] = fit_naked_dna(Turn2_mean_1way,Z_beadc2_mean_1way(:,j));

            % Height at (+) and (-) buckling 
            slope_left = abs(paraf(7) ); % um/turn
            slope_right = abs(paraf(6) ); % um/turn
            trace{j}.slope_left = slope_left;
            trace{j}.slope_right = slope_right;
            trace{j}.height_0 = height_0;

            %% Analyze topo activity, deduce wait time and average rate
            %[sign_winding1, t_on_topo1, z_holding1_smooth,index_hold1, index_no_topo_hold1, v_relaxed1, index_no_topo1,index_topo_activity_end_excluded1] = ...
            %    analayze_topo_activity(time_holding1,turn_holding1,z_holding1 , paraf, fitresult_height_turn_right, fitresult_height_turn_left);

            %[sign_winding2, t_on_topo2, z_holding2_smooth, index_hold2, index_no_topo_hold2, v_relaxed2, index_no_topo2,index_topo_activity_end_excluded2] = ...
            %    analayze_topo_activity(time_holding2,turn_holding2,z_holding2 , paraf, fitresult_height_turn_right, fitresult_height_turn_left);

            %[sign_winding3, t_on_topo3, z_holding3_smooth, index_hold3, index_no_topo_hold3, v_relaxed3, index_no_topo3,index_topo_activity_end_excluded3] = ...
            %    analayze_topo_activity(time_holding3,turn_holding3,z_holding3 , paraf, fitresult_height_turn_right, fitresult_height_turn_left);

            [Z_holding_s1_1, ~] = ...
            analayze_topo_activity(time_holding1,turn_holding1,z_holding1, trace{j}.height_c_before);
            [Z_holding_s2_1, ~] = ...
            analayze_topo_activity(time_holding2,turn_holding2,z_holding2);
            %[v_pausefree_3, dwelltime_pause_3, sign_winding_3, t_on_topo_3, Z_holding_s1_3, Z_holding_s2_3, v_3, time_3, Z_3, turn_3, index_pause_time_3, pause_density_3, distance_between_pauses_3, v_pausefree_mean_3, index_hold_3, index_no_topo_hold_3] = ...
            %    analayze_topo_activity(time_holding3,turn_holding3,z_holding3);
            
            %trace{j}.t_on_topo = [t_on_topo_1 t_on_topo_2 t_on_topo_3 ];
%             trace{j}.t_on_topo = t_on_topo_1 ;
%             %trace{j}.sign_winding =  [sign_winding_1 sign_winding_2 sign_winding_3 ];
%             trace{j}.sign_winding =  sign_winding_1 ;
%             %trace{j}.dwelltime_pause =  {dwelltime_pause_1 dwelltime_pause_2 dwelltime_pause_3 };
%             trace{j}.dwelltime_pause =  dwelltime_pause_1 ;
%             
%             trace{j}.v_mean_1 = v_mean_1;
%             trace{j}.Z_pausefree_1 = Z_pausefree_1;
%             trace{j}.Z_topo_pause_mean_1 = Z_topo_pause_mean_1;
%             
%             t_pause_threshold = 10;
%             num_pause_buckling = sum(Z_topo_pause_mean_1 <=  trace{j}.height_c_before * 0.9  & dwelltime_pause_1 >= t_pause_threshold);
%             trace{j}.num_pause_buckling = num_pause_buckling;
%             
%             %trace{j}.distance_between_pauses =  {distance_between_pauses_1 distance_between_pauses_2 distance_between_pauses_3 };
%             trace{j}.distance_between_pauses =  distance_between_pauses_1 ;
%             
%             distance_between_pauses_turn = distance_between_pauses_1 / slope_left;
%             trace{j}.distance_between_pauses_turn =  distance_between_pauses_turn;
%             
%             %trace{j}.v_pausefree_mean =  {v_pausefree_mean_1 v_pausefree_mean_2 v_pausefree_mean_3 };
%             trace{j}.v_pausefree_mean =  v_pausefree_mean_1 ;
%             
%             v_pausefree_mean_tps = v_pausefree_mean_1/slope_left;
%             trace{j}.v_pausefree_mean_tps =  v_pausefree_mean_tps;
%             
%             
%             %trace{j}.v_pausefree =  {v_pausefree_1 v_pausefree_2 v_pausefree_3 };
%             trace{j}.v_pausefree =  v_pausefree_1;
%             v_pausefree_tps = v_pausefree_1/slope_left;
%             trace{j}.v_pausefree_tps =  v_pausefree_tps;
%             
%             %trace{j}.pause_density =  {pause_density_1 pause_density_2 pause_density_3 };
%             trace{j}.pause_density =  pause_density_1;
%             
            %% plot data
            f2 = figure;
            pos = [10 70 900 950];
            set(f2, 'Pos', pos);
            subplot(4,1,1:2);
            hold all
            
            line([0 0],[ -2 15],'LineStyle','--','LineWidth',1, 'Color','k');
            line([d_wound1 d_wound1],[ -2 15],'LineStyle','--','LineWidth',1, 'Color',[0.87, 0.49, 0]);
            %line([d_wound2 d_wound2],[ -2 15],'LineStyle','--','LineWidth',1, 'Color','b');
            %line([d_wound3 d_wound3],[ -2 15],'LineStyle','--','LineWidth',1, 'Color',[0, 0.5, 0]);
            plot( ttempt, ztempt,'-c','LineWidth',2);
            
            plot(Turn4_mean_1way, Z_beadc4_mean_1way(:,j), '-o','MarkerSize',2,'LineWidth',1,'MarkerFaceColor',[0.5, 0.5, 0.5],'Color',[0.5, 0.5, 0.5]);
            %plot(Turn4_mean_1way_negative, Z_beadc4_mean_1way_negative(:,j), '-o','MarkerSize',2,'LineWidth',1,'MarkerFaceColor',[0.5, 0.5, 0.5],'Color',[0.5, 0.5, 0.5]);
        
            plot(movingmean(Turn3,10), movingmean(Z_beadc3(:,j),10), '.','Color',[0.87, 0.49, 0]);
            %plot(movingmean(Turn5,10), movingmean(Z_beadc5(:,j),10), '.','Color','b');
            %plot(movingmean(Turn6,10), movingmean(Z_beadc6(:,j),10), '.','Color',[0, 0.5, 0]);
            
            
            plot(Turn2_mean_1way, Z_beadc2_mean_1way(:,j), '-ok','MarkerSize',4,'LineWidth',1.5,'MarkerFaceColor','k');
            %plot(Turn7_mean_1way, Z_beadc7_mean_1way(:,j), '-o','MarkerSize',4,'LineWidth',1.5,'MarkerFaceColor','r','Color','r');
            

            force_topo = num2str(F_from_var_length_p);
            title([Box_namec1{j} ', F = '  force_topo(1:4) 'pN'],'Interpreter','None','FontSize',12,'FontWeight','bold');
            ylabel('Z (\mum)','FontSize',12);
            xlabel('Turn added','FontSize',12);
            set(gca,'FontSize',12,'FontName','Calibri');
            hmm = 1.1*max(Z_beadc2_mean_1way(:,j));

            xlim([-60 150]); ylim([0 hmm+0.2]);
            textp1 = { ['Max height: ' , num2str(height_0) ' \mum'] };   
            %text(Turn4(end)+4,hmm+0.19, textp1, 'HorizontalAlignment', 'right','VerticalAlignment', 'top','FontSize',12, 'FontWeight','bold', 'BackgroundColor','w', 'EdgeColor','k');
            
            
            %t_no_topo = time_holding1(index_no_topo_hold_1);
            
            
            subplot(4,1,3);
            hold all
            line([0 time_holding1(end)+5],[height_0 height_0],'LineStyle','--','LineWidth',1, 'Color',[0.5, 0.5, 0.5]);
            plot(time_holding1, Z_holding_s1_1, ':','Color',[0, 0.75, 0.75], 'LineWidth',1);
            %plot(time_holding1(index_hold_1), Z_holding_s1_1(index_hold_1), '-','Color',[0, 0.75, 0.75], 'LineWidth',1);
            
            %plot(time_pausefree_squeeze_1 + t_no_topo(end), Z_pausefree_1, '.--' ,'Color','r', 'LineWidth',1, 'MarkerSize',3);
            
            %plot(time_1(cell2mat(index_pause_time_1(:))), Z_1(cell2mat(index_pause_time_1(:))),'.k', 'MarkerSize',2);
            %plot(time_holding1(index_no_topo_hold_1), Z_holding_s1_1(index_no_topo_hold_1),'-','LineWidth',0.5, 'Color',[0.5,0.5,0.5]);
            xlabel('Time (s)', 'FontSize',12); 
            ylabel('Height (\mum)', 'FontSize',12);
            xlim([0 920]); ylim([min(Z_holding_s1_1)-0.1 hmm+0.1]);
            %title('Positive side','FontSize',13);
            set(gca,'FontSize',12,'FontName','Calibri');     
            
            %title(['On-time: ' num2str(t_on_topo_1) ' s, pause-free rate: ' num2str(v_pausefree_mean_1) ' um/s, ' num2str(v_pausefree_mean_tps) ' tps, # of long pauses after buckling: ' num2str(num_pause_buckling)],'Interpreter','None','FontSize',12,'FontWeight','bold');

            subplot(4,1,4);
            hold all
            line([0 time_holding2(end)+5],[height_0 height_0],'LineStyle','--','LineWidth',1, 'Color',[0.5, 0.5, 0.5]);
            plot(time_holding2, Z_holding_s2_1, '-','Color','b', 'LineWidth',1.5);

%             if any(index_no_topo_hold2)
%                 plot(time_holding2(index_no_topo_hold2), z_holding2_smooth(index_no_topo_hold2), '+','Color','k');
%             end
% 
%             if any(index_topo_activity_end_excluded2)
%                 plot(time_holding2(index_topo_activity_end_excluded2), z_holding2_smooth(index_topo_activity_end_excluded2), '+','Color','r');
%             end
            ylabel('Height (\mum)', 'FontSize',12);
            xlim([0 1820]); ylim([0 hmm+0.2]);
            set(gca,'FontSize',12,'FontName','Calibri');     
            
            
            prompt = {'Comment?'};
            topheader = 'Add comment';
            num_line = 1;
            defaultvals = {'good'};
            comments = inputdlg(prompt, topheader,num_line,defaultvals);
            trace{j}.comments = comments;        

            savefig(f2,[Box_namec1{j} '.fig']);
            %print(Box_namec1{j},'-dpng','-r0');
            %export_fig(f2,[Box_namec1{j} '.png']);


            %data1 = [trace{j}.F_from_var_length_p height_0  sign_winding_1 t_on_topo_1 sign_winding_2 t_on_topo_2 sign_winding_3 t_on_topo_3 ];
            data1 = [trace{j}.F_from_var_length_p height_0] % sign_winding_1 t_on_topo_1     ];
            data1 =  [trace{j}.name trace{j}.comments num2cell(data1)];
            if ~isempty(data1)
                        xlswrite('summary.xls', data1, 1, ['A' num2str(j)]);
            end

        else
            prompt = {'Comment?'};
            topheader = 'Add comment';
            num_line = 1;
            defaultvals = {'bad_hysteresis'};
            comments = inputdlg(prompt, topheader,num_line,defaultvals);
            trace{j}.comments = comments;
            data1 = [trace{j}.name comments];
            if ~isempty(data1)
                        xlswrite('summary.xls', data1, 2, ['A' num2str(j)]);
            end
            close(ftemp);
        end
    end
    
    save('data.mat','trace');
end
    
%save('data.mat','trace');

end

function [Turn2_mean_1way, Z_beadc2_mean_1way, Turn2_mean_1way_negative, Z_beadc2_mean_1way_negative] = average_at_pause(Times2, Turn2, Z_beadc2)
    k=1;
    time_pat{1} = Times2(1);
    turn_pat{1} = Turn2(1);
    index_pat{1} = 1;
    for i = 1 : length(Turn2)-1
            if abs(Turn2(i+1)-turn_pat{k}(end))<0.01
                time_pat{k} = [time_pat{k} Times2(i+1)];   
                turn_pat{k} = [turn_pat{k} Turn2(i+1)];   
                index_pat{k} = [index_pat{k} i+1];   
            else
                if abs(Turn2(i+1) - turn_pat{k}(end)) == 1
                    k = k+1;
                    time_pat{k} = Times2(i+1);
                    turn_pat{k} = Turn2(i+1);
                    index_pat{k} = i+1;
                end
            end
    end    
    n1 = length(time_pat);
    Turn2_mean  = zeros(1,n1);
    L = size(Z_beadc2,2);
    Z_beadc2_mean  = zeros(n1,L);
    for i = 1 : n1
            Turn2_mean(i) = mean(turn_pat{i}); 
            Z_beadc2_mean(i,:) = mean(Z_beadc2(index_pat{i},:),1);
    end
    index_pos_step2 = diff(Turn2_mean)>0; index_pos_step2_negative = diff(Turn2_mean)<0;
    Turn2_mean_1way = Turn2_mean(index_pos_step2); Z_beadc2_mean_1way = Z_beadc2_mean(index_pos_step2,:);
    Turn2_mean_1way_negative = Turn2_mean(index_pos_step2_negative); Z_beadc2_mean_1way_negative = Z_beadc2_mean(index_pos_step2_negative,:);
end


function [height_c_before, F_from_var_length, F_from_var_length_p] = force_calculation(Times1, X_beadc1, Y_beadc1, Z_beadc1, uB, height_correction, tv)
        
        kbt = 4.08;
        R_bead = 0.5;
        height_c_before = mean(Z_beadc1);
               
        index = find(Times1<=tv);
        [xData, yData] = prepareCurveData( Times1(index), X_beadc1(index));
        
        ft1 = fittype( 'poly1' );
        [fitresult1, ~] = fit( xData, yData, ft1 );
        xx = X_beadc1(index) - fitresult1(Times1(index));
        [xData, yData] = prepareCurveData( Times1(index), Y_beadc1(index));
        ft2 = fittype( 'poly1' );
        [fitresult2, ~] = fit( xData, yData, ft2 );
        yy = Y_beadc1(index) - fitresult2(Times1(index));
        uBp = [-uB(2), uB(1)]; uBp = uBp/norm(uBp);
        rB = [xx yy] * uB'; % projection on the B field
        rBp = [xx yy] * uBp'; % projection on the perpendicular direction to B field
        uB_var1 = mean(rB.^2) - mean(rB)^2;    % variance along along B   
        uB_var2 = mean(rBp.^2) - mean(rBp)^2;

        F_from_var_length = kbt * height_c_before*1e3/(uB_var1*1e6);                        % parallel force
        F_from_var_length_p = kbt * (height_c_before + R_bead - height_correction)*1e3/(uB_var2*1e6);    % perpendicular force 
end

function [paraf, height_0, turn_0, ttempt, ztempt,fitresult_height_turn_right, fitresult_height_turn_left] = fit_naked_dna(turn_c,z_m_c)
        % use 3-piece function to fit the data for dna % get the residue of fitting
        range_f = [-60 60];
        turn_c = turn_c(:);
        z_m_c = z_m_c(:);
        turn_c_f = turn_c(turn_c >= range_f(1) & turn_c <= range_f(2));
        z_m_c_f = z_m_c(turn_c >= range_f(1) & turn_c <= range_f(2));
        
        index_max = find(z_m_c_f >= max(z_m_c_f)*0.9);
        xf = turn_c_f( index_max);
        yf = z_m_c_f( index_max);
        [xData, yData] = prepareCurveData( xf, yf );
        ft = fittype( 'poly2' );
        [fitparabol, ~] = fit( xData, yData, ft );
        turn_guess = -fitparabol.p2 / (2 * fitparabol.p1);
        h_max_guess = fitparabol.p3 - fitparabol.p2^2/(4*fitparabol.p1);

        % for DNA, piecewise
        [xData, yData] = prepareCurveData( turn_c_f, z_m_c_f );       
        F_positive = @(para, x) (para(1) + para(2) * (x - para(3)).^2) .* (x < para(4) & x > para(5)) + (para(1) + para(2)* (para(4) - para(3)).^2 +para(6).*(x-para(4))) .* (x > para(4)) + ...
                            + (para(1) + para(2)* (para(5) - para(3)).^2 +para(7).*(x - para(5))) .* (x < para(5)) ;
        para0 = [h_max_guess, -0.001 , turn_guess, turn_guess+13, turn_guess-13, -60/1000, 60/1000];
        [paraf,~,~,~,~, ~, ~] = lsqcurvefit(F_positive,para0,xData, yData);
        
        height_0 = paraf(1);
        
        %height_0 = mean(z_m_c_f(abs(turn_c_f - turn_c_f(z_m_c_f == max(z_m_c_f))) <= 2));
        turn_0 = paraf(3);
        
        ttempt = -65:0.1:65;
        ztempt = F_positive(paraf,ttempt);
        
        % positive turns
        fitresult_height_turn_right = @(para, h) (para(4) + 1/para(6) * (h - para(1) - para(2) * (para(4)- para(3))^2) ) .* (h < para(1)+para(2)*(para(4) - para(3))^2) + ...
                    (para(3) + sqrt( (h - para(1))/para(2) ) ) .* (h >= para(1)+para(2)*(para(4)-para(3))^2 & h < para(1)) + para(3)*( h >= para(1));
        % negative turns
        fitresult_height_turn_left = @(para, h) (para(5) + 1/para(7) * (h - para(1) - para(2) * (para(5)- para(3))^2) ) .* (h < para(1)+para(2)*(para(5) - para(3))^2) + ...
                    (para(3) - sqrt( (h - para(1))/para(2) ) ) .* (h >= para(1)+para(2)*(para(5)-para(3))^2 & h < para(1)) + para(3)*( h >= para(1));
end

function [Z_holding_s1, Z_holding_s2] = analayze_topo_activity(time_holding,turn_holding, Z_holding, h_max)  
        d_wound = turn_holding(end);
        d_0 = turn_holding(1);
        sign_winding = mean(diff(turn_holding))>0; % positive if winding positively
        
        rate = 1/0.4;  
        dt = mean(diff(time_holding));
        window = round(1/rate / dt);
        order = 2;
        if mod(window, 2) == 0
                window = window+1;
        end
        Z_holding_s1 = sgolayfilt(Z_holding,order,window);
         
        window2 = round(2.5/ dt);
        if mod(window2, 2) == 0
                window2 = window2+1;
        end
        Z_holding_s2 = sgolayfilt(Z_holding,order,window2);
end

function Z_beadc3 = z_correct(Times3, Z_beadc3, L)

zmin = -1; zmax = 6;
for i = 1 :L
    ztemp = Z_beadc3(:,i);
    index_large = ztemp<zmin | ztemp>zmax;
    ttemp = Times3(~index_large);
    ztemptt = ztemp(~index_large);
    
    if sum(index_large) <= size(ttemp,1)
        if mean(ztemp)> zmin && mean(ztemp)< zmax
            if sum(index_large) >=1
                ztemp(index_large) = interp1(ttemp(:), ztemptt(:), Times3(index_large));
                Z_beadc3(:,i) = ztemp;
            end
        end
    end
end
end

function [Times, Turn, Box_name, X_bead, Y_bead, Z_bead,L, Z_piezo] = get_mtdata4(filename)
% use raw data without correction

A = importdata(filename);
Times = A.data(:,1)/1000;           %s 
    Turn = A.data(:, 2)/32; Turn_c = A.data(:, 7)/32; Turn = Turn + Turn_c(1); %turn
    Box_data = A.data(:,15:end);
    Box_name = A.colheaders(15:3:end);
    L = length(Box_name);
    X_bead = Box_data(:,1:3:end) * 0.1465 * 2;
    Y_bead = Box_data(:,2:3:end) * 0.1465 * 2;
    Z_bead = Box_data(:,3:3:end) * 1.2*1.33;
    Z_piezo = A.data(:,4)  * 1.2*1.33;    
end