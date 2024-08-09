% Plotting PSTH in different saccades

load("SNr_new.mat")

% Assuming dataTable is your table name and SaccNUM is the column of interest
uniqueValues = unique(table.SaccNum);
numberOfUniqueValues = length(uniqueValues);

% Display the number of unique values
disp(['Number of unique values in SaccNUM: ', num2str(numberOfUniqueValues)]);
% The number of different Sacc numbers are 11

column_names = cell(1, 1600);

for i = 1:1600
    column_names{i} = ['bin', num2str(i)];
end




%% Arranging Data

% TP/TA separation
TempTP=find(table.EventValue==4);
TempTA=find(table.EventValue==3);

% Eff/Ineff separation
Eff=find(table.Search_Type==1);
Ineff=find(table.Search_Type==0);

%Efficient=find(tableB.SlopeTP<=20);
%Inefficient=find(tableB.SlopeTP>=35);


% Sacc Numbers separation
Sacc1=find(table.SaccNum==1);
Sacc2=find(table.SaccNum==2);
Sacc3=find(table.SaccNum==3);
Sacc4=find(table.SaccNum==4);
Sacc5=find(table.SaccNum==5);
Sacc6=find(table.SaccNum==6);
Sacc7=find(table.SaccNum==7);
Sacc8=find(table.SaccNum==8);
Sacc9=find(table.SaccNum==9);
Sacc10=find(table.SaccNum==10);
Sacc11=find(table.SaccNum==11);


% 1 Saccade
TP_ALL_Eff_Sacc1_temp=intersect(intersect(TempTP,Eff),Sacc1);
TP_ALL_Ineff_Sacc1_temp=intersect(intersect(TempTP,Ineff),Sacc1);
TA_ALL_Eff_Sacc1_temp=intersect(intersect(TempTA,Eff),Sacc1);
TA_ALL_Ineff_Sacc1_temp=intersect(intersect(TempTA,Ineff),Sacc1);
TA_ALL_Ineff_Sacc1=table(TA_ALL_Ineff_Sacc1_temp,:);
TP_ALL_Ineff_Sacc1=table(TP_ALL_Ineff_Sacc1_temp,:);
TP_ALL_Eff_Sacc1=table(TP_ALL_Eff_Sacc1_temp,:);
TA_ALL_Eff_Sacc1=table(TA_ALL_Eff_Sacc1_temp,:);
TA_ALL_Eff_Sacc1_Bins=TA_ALL_Eff_Sacc1{:, column_names(1:1600)};
TA_ALL_Ineff_Sacc1_Bins=TA_ALL_Ineff_Sacc1{:, column_names(1:1600)};
TP_ALL_Ineff_Sacc1_Bins=TP_ALL_Ineff_Sacc1{:, column_names(1:1600)};
TP_ALL_Eff_Sacc1_Bins=TP_ALL_Eff_Sacc1{:, column_names(1:1600)};
TP_ALL_Eff_Sacc1_Avg=nanmean(TP_ALL_Eff_Sacc1_Bins,1);
TP_ALL_Ineff_Sacc1_Avg=nanmean(TP_ALL_Ineff_Sacc1_Bins,1);
TA_ALL_Ineff_Sacc1_Avg=nanmean(TA_ALL_Ineff_Sacc1_Bins,1);
TA_ALL_Eff_Sacc1_Avg=nanmean(TA_ALL_Eff_Sacc1_Bins,1);

% 2 Saccade
TP_ALL_Eff_Sacc2_temp=intersect(intersect(TempTP,Eff),Sacc2);
TP_ALL_Ineff_Sacc2_temp=intersect(intersect(TempTP,Ineff),Sacc2);
TA_ALL_Eff_Sacc2_temp=intersect(intersect(TempTA,Eff),Sacc2);
TA_ALL_Ineff_Sacc2_temp=intersect(intersect(TempTA,Ineff),Sacc2);
TA_ALL_Ineff_Sacc2=table(TA_ALL_Ineff_Sacc2_temp,:);
TP_ALL_Ineff_Sacc2=table(TP_ALL_Ineff_Sacc2_temp,:);
TP_ALL_Eff_Sacc2=table(TP_ALL_Eff_Sacc2_temp,:);
TA_ALL_Eff_Sacc2=table(TA_ALL_Eff_Sacc2_temp,:);
TA_ALL_Eff_Sacc2_Bins=TA_ALL_Eff_Sacc2{:, column_names(1:1600)};
TA_ALL_Ineff_Sacc2_Bins=TA_ALL_Ineff_Sacc2{:, column_names(1:1600)};
TP_ALL_Ineff_Sacc2_Bins=TP_ALL_Ineff_Sacc2{:, column_names(1:1600)};
TP_ALL_Eff_Sacc2_Bins=TP_ALL_Eff_Sacc2{:, column_names(1:1600)};
TP_ALL_Eff_Sacc2_Avg=nanmean(TP_ALL_Eff_Sacc2_Bins,1);
TP_ALL_Ineff_Sacc2_Avg=nanmean(TP_ALL_Ineff_Sacc2_Bins,1);
TA_ALL_Ineff_Sacc2_Avg=nanmean(TA_ALL_Ineff_Sacc2_Bins,1);
TA_ALL_Eff_Sacc2_Avg=nanmean(TA_ALL_Eff_Sacc2_Bins,1);

% 3 Saccade
TP_ALL_Eff_Sacc3_temp=intersect(intersect(TempTP,Eff),Sacc3);
TP_ALL_Ineff_Sacc3_temp=intersect(intersect(TempTP,Ineff),Sacc3);
TA_ALL_Eff_Sacc3_temp=intersect(intersect(TempTA,Eff),Sacc3);
TA_ALL_Ineff_Sacc3_temp=intersect(intersect(TempTA,Ineff),Sacc3);
TA_ALL_Ineff_Sacc3=table(TA_ALL_Ineff_Sacc3_temp,:);
TP_ALL_Ineff_Sacc3=table(TP_ALL_Ineff_Sacc3_temp,:);
TP_ALL_Eff_Sacc3=table(TP_ALL_Eff_Sacc3_temp,:);
TA_ALL_Eff_Sacc3=table(TA_ALL_Eff_Sacc3_temp,:);
TA_ALL_Eff_Sacc3_Bins=TA_ALL_Eff_Sacc3{:, column_names(1:1600)};
TA_ALL_Ineff_Sacc3_Bins=TA_ALL_Ineff_Sacc3{:, column_names(1:1600)};
TP_ALL_Ineff_Sacc3_Bins=TP_ALL_Ineff_Sacc3{:, column_names(1:1600)};
TP_ALL_Eff_Sacc3_Bins=TP_ALL_Eff_Sacc3{:, column_names(1:1600)};
TP_ALL_Eff_Sacc3_Avg=nanmean(TP_ALL_Eff_Sacc3_Bins,1);
TP_ALL_Ineff_Sacc3_Avg=nanmean(TP_ALL_Ineff_Sacc3_Bins,1);
TA_ALL_Ineff_Sacc3_Avg=nanmean(TA_ALL_Ineff_Sacc3_Bins,1);
TA_ALL_Eff_Sacc3_Avg=nanmean(TA_ALL_Eff_Sacc3_Bins,1);

% 4 Saccade
TP_ALL_Eff_Sacc4_temp=intersect(intersect(TempTP,Eff),Sacc4);
TP_ALL_Ineff_Sacc4_temp=intersect(intersect(TempTP,Ineff),Sacc4);
TA_ALL_Eff_Sacc4_temp=intersect(intersect(TempTA,Eff),Sacc4);
TA_ALL_Ineff_Sacc4_temp=intersect(intersect(TempTA,Ineff),Sacc4);
TA_ALL_Ineff_Sacc4=table(TA_ALL_Ineff_Sacc4_temp,:);
TP_ALL_Ineff_Sacc4=table(TP_ALL_Ineff_Sacc4_temp,:);
TP_ALL_Eff_Sacc4=table(TP_ALL_Eff_Sacc4_temp,:);
TA_ALL_Eff_Sacc4=table(TA_ALL_Eff_Sacc4_temp,:);
TA_ALL_Eff_Sacc4_Bins=TA_ALL_Eff_Sacc4{:, column_names(1:1600)};
TA_ALL_Ineff_Sacc4_Bins=TA_ALL_Ineff_Sacc4{:, column_names(1:1600)};
TP_ALL_Ineff_Sacc4_Bins=TP_ALL_Ineff_Sacc4{:, column_names(1:1600)};
TP_ALL_Eff_Sacc4_Bins=TP_ALL_Eff_Sacc4{:, column_names(1:1600)};
TP_ALL_Eff_Sacc4_Avg=nanmean(TP_ALL_Eff_Sacc4_Bins,1);
TP_ALL_Ineff_Sacc4_Avg=nanmean(TP_ALL_Ineff_Sacc4_Bins,1);
TA_ALL_Ineff_Sacc4_Avg=nanmean(TA_ALL_Ineff_Sacc4_Bins,1);
TA_ALL_Eff_Sacc4_Avg=nanmean(TA_ALL_Eff_Sacc4_Bins,1);


% 5 Saccade
TP_ALL_Eff_Sacc5_temp=intersect(intersect(TempTP,Eff),Sacc5);
TP_ALL_Ineff_Sacc5_temp=intersect(intersect(TempTP,Ineff),Sacc5);
TA_ALL_Eff_Sacc5_temp=intersect(intersect(TempTA,Eff),Sacc5);
TA_ALL_Ineff_Sacc5_temp=intersect(intersect(TempTA,Ineff),Sacc5);
TA_ALL_Ineff_Sacc5=table(TA_ALL_Ineff_Sacc5_temp,:);
TP_ALL_Ineff_Sacc5=table(TP_ALL_Ineff_Sacc5_temp,:);
TP_ALL_Eff_Sacc5=table(TP_ALL_Eff_Sacc5_temp,:);
TA_ALL_Eff_Sacc5=table(TA_ALL_Eff_Sacc5_temp,:);
TA_ALL_Eff_Sacc5_Bins=TA_ALL_Eff_Sacc5{:, column_names(1:1600)};
TA_ALL_Ineff_Sacc5_Bins=TA_ALL_Ineff_Sacc5{:, column_names(1:1600)};
TP_ALL_Ineff_Sacc5_Bins=TP_ALL_Ineff_Sacc5{:, column_names(1:1600)};
TP_ALL_Eff_Sacc5_Bins=TP_ALL_Eff_Sacc5{:, column_names(1:1600)};
TP_ALL_Eff_Sacc5_Avg=nanmean(TP_ALL_Eff_Sacc5_Bins,1);
TP_ALL_Ineff_Sacc5_Avg=nanmean(TP_ALL_Ineff_Sacc5_Bins,1);
TA_ALL_Ineff_Sacc5_Avg=nanmean(TA_ALL_Ineff_Sacc5_Bins,1);
TA_ALL_Eff_Sacc5_Avg=nanmean(TA_ALL_Eff_Sacc5_Bins,1);

% 6 Saccade
TP_ALL_Eff_Sacc6_temp=intersect(intersect(TempTP,Eff),Sacc6);
TP_ALL_Ineff_Sacc6_temp=intersect(intersect(TempTP,Ineff),Sacc6);
TA_ALL_Eff_Sacc6_temp=intersect(intersect(TempTA,Eff),Sacc6);
TA_ALL_Ineff_Sacc6_temp=intersect(intersect(TempTA,Ineff),Sacc6);
TA_ALL_Ineff_Sacc6=table(TA_ALL_Ineff_Sacc6_temp,:);
TP_ALL_Ineff_Sacc6=table(TP_ALL_Ineff_Sacc6_temp,:);
TP_ALL_Eff_Sacc6=table(TP_ALL_Eff_Sacc6_temp,:);
TA_ALL_Eff_Sacc6=table(TA_ALL_Eff_Sacc6_temp,:);
TA_ALL_Eff_Sacc6_Bins=TA_ALL_Eff_Sacc6{:, column_names(1:1600)};
TA_ALL_Ineff_Sacc6_Bins=TA_ALL_Ineff_Sacc6{:, column_names(1:1600)};
TP_ALL_Ineff_Sacc6_Bins=TP_ALL_Ineff_Sacc6{:, column_names(1:1600)};
TP_ALL_Eff_Sacc6_Bins=TP_ALL_Eff_Sacc6{:, column_names(1:1600)};
TP_ALL_Eff_Sacc6_Avg=nanmean(TP_ALL_Eff_Sacc6_Bins,1);
TP_ALL_Ineff_Sacc6_Avg=nanmean(TP_ALL_Ineff_Sacc6_Bins,1);
TA_ALL_Ineff_Sacc6_Avg=nanmean(TA_ALL_Ineff_Sacc6_Bins,1);
TA_ALL_Eff_Sacc6_Avg=nanmean(TA_ALL_Eff_Sacc6_Bins,1);

% 7 Saccade
TP_ALL_Eff_Sacc7_temp=intersect(intersect(TempTP,Eff),Sacc7);
TP_ALL_Ineff_Sacc7_temp=intersect(intersect(TempTP,Ineff),Sacc7);
TA_ALL_Eff_Sacc7_temp=intersect(intersect(TempTA,Eff),Sacc7);
TA_ALL_Ineff_Sacc7_temp=intersect(intersect(TempTA,Ineff),Sacc7);
TA_ALL_Ineff_Sacc7=table(TA_ALL_Ineff_Sacc7_temp,:);
TP_ALL_Ineff_Sacc7=table(TP_ALL_Ineff_Sacc7_temp,:);
TP_ALL_Eff_Sacc7=table(TP_ALL_Eff_Sacc7_temp,:);
TA_ALL_Eff_Sacc7=table(TA_ALL_Eff_Sacc7_temp,:);
TA_ALL_Eff_Sacc7_Bins=TA_ALL_Eff_Sacc7{:, column_names(1:1600)};
TA_ALL_Ineff_Sacc7_Bins=TA_ALL_Ineff_Sacc7{:, column_names(1:1600)};
TP_ALL_Ineff_Sacc7_Bins=TP_ALL_Ineff_Sacc7{:, column_names(1:1600)};
TP_ALL_Eff_Sacc7_Bins=TP_ALL_Eff_Sacc7{:, column_names(1:1600)};
TP_ALL_Eff_Sacc7_Avg=nanmean(TP_ALL_Eff_Sacc7_Bins,1);
TP_ALL_Ineff_Sacc7_Avg=nanmean(TP_ALL_Ineff_Sacc7_Bins,1);
TA_ALL_Ineff_Sacc7_Avg=nanmean(TA_ALL_Ineff_Sacc7_Bins,1);
TA_ALL_Eff_Sacc7_Avg=nanmean(TA_ALL_Eff_Sacc7_Bins,1);

% 8 Saccade
TP_ALL_Eff_Sacc8_temp=intersect(intersect(TempTP,Eff),Sacc8);
TP_ALL_Ineff_Sacc8_temp=intersect(intersect(TempTP,Ineff),Sacc8);
TA_ALL_Eff_Sacc8_temp=intersect(intersect(TempTA,Eff),Sacc8);
TA_ALL_Ineff_Sacc8_temp=intersect(intersect(TempTA,Ineff),Sacc8);
TA_ALL_Ineff_Sacc8=table(TA_ALL_Ineff_Sacc8_temp,:);
TP_ALL_Ineff_Sacc8=table(TP_ALL_Ineff_Sacc8_temp,:);
TP_ALL_Eff_Sacc8=table(TP_ALL_Eff_Sacc8_temp,:);
TA_ALL_Eff_Sacc8=table(TA_ALL_Eff_Sacc8_temp,:);
TA_ALL_Eff_Sacc8_Bins=TA_ALL_Eff_Sacc8{:, column_names(1:1600)};
TA_ALL_Ineff_Sacc8_Bins=TA_ALL_Ineff_Sacc8{:, column_names(1:1600)};
TP_ALL_Ineff_Sacc8_Bins=TP_ALL_Ineff_Sacc8{:, column_names(1:1600)};
TP_ALL_Eff_Sacc8_Bins=TP_ALL_Eff_Sacc8{:, column_names(1:1600)};
TP_ALL_Eff_Sacc8_Avg=nanmean(TP_ALL_Eff_Sacc8_Bins,1);
TP_ALL_Ineff_Sacc8_Avg=nanmean(TP_ALL_Ineff_Sacc8_Bins,1);
TA_ALL_Ineff_Sacc8_Avg=nanmean(TA_ALL_Ineff_Sacc8_Bins,1);
TA_ALL_Eff_Sacc8_Avg=nanmean(TA_ALL_Eff_Sacc8_Bins,1);


%% Plotting for Eff

time_axis = linspace(-0.6, 1.0, 1600);
colorsC = get_distinguishable_colors(2); % You can also use other colormap functions like parula, jet, etc.
 

figure; % Creates a new figure window
 % Create subplot in a 2x4 grid at the i-th position
 %sacc1
subplot(2, 4, 1);
hold on;
PlotPSTH(time_axis,TP_ALL_Eff_Sacc1_Avg,colorsC(1,:));
PlotPSTH(time_axis,TA_ALL_Eff_Sacc1_Avg,colorsC(2,:));
% Labeling
xlabel('Time (s)');
ylabel('Firing Rate (spikes/s)');
title('Efficient search wirh SaccNum=1');
legend('TP', 'TA');
grid on;
hold off;

%sacc2
subplot(2, 4, 2);
hold on;
PlotPSTH(time_axis,TP_ALL_Eff_Sacc2_Avg,colorsC(1,:));
PlotPSTH(time_axis,TA_ALL_Eff_Sacc2_Avg,colorsC(2,:));
% Labeling
xlabel('Time (s)');
ylabel('Firing Rate (spikes/s)');
title('Efficient search wirh SaccNum=2');
legend('TP', 'TA');
grid on;
hold off;

%sacc3
subplot(2, 4, 3);
hold on;
PlotPSTH(time_axis,TP_ALL_Eff_Sacc3_Avg,colorsC(1,:));
PlotPSTH(time_axis,TA_ALL_Eff_Sacc3_Avg,colorsC(2,:));
% Labeling
xlabel('Time (s)');
ylabel('Firing Rate (spikes/s)');
title('Efficient search wirh SaccNum=3');
legend('TP', 'TA');
grid on;
hold off;

%sacc4
subplot(2, 4, 4);
hold on;
PlotPSTH(time_axis,TP_ALL_Eff_Sacc4_Avg,colorsC(1,:));
PlotPSTH(time_axis,TA_ALL_Eff_Sacc4_Avg,colorsC(2,:));
% Labeling
xlabel('Time (s)');
ylabel('Firing Rate (spikes/s)');
title('Efficient search wirh SaccNum=4');
legend('TP', 'TA');
grid on;
hold off;

%sacc5
subplot(2, 4, 5);
hold on;
PlotPSTH(time_axis,TP_ALL_Eff_Sacc5_Avg,colorsC(1,:));
PlotPSTH(time_axis,TA_ALL_Eff_Sacc5_Avg,colorsC(2,:));
% Labeling
xlabel('Time (s)');
ylabel('Firing Rate (spikes/s)');
title('Efficient search wirh SaccNum=5');
legend('TP', 'TA');
grid on;
hold off;

%sacc6
subplot(2, 4, 6);
hold on;
PlotPSTH(time_axis,TP_ALL_Eff_Sacc6_Avg,colorsC(1,:));
PlotPSTH(time_axis,TA_ALL_Eff_Sacc6_Avg,colorsC(2,:));
% Labeling
xlabel('Time (s)');
ylabel('Firing Rate (spikes/s)');
title('Efficient search wirh SaccNum=6');
legend('TP', 'TA');
grid on;
hold off;

%sacc7
subplot(2, 4, 7);
hold on;
PlotPSTH(time_axis,TP_ALL_Eff_Sacc7_Avg,colorsC(1,:));
PlotPSTH(time_axis,TA_ALL_Eff_Sacc7_Avg,colorsC(2,:));
% Labeling
xlabel('Time (s)');
ylabel('Firing Rate (spikes/s)');
title('Efficient search wirh SaccNum=7');
legend('TP', 'TA');
grid on;
hold off;

%sacc8
subplot(2, 4, 8);
hold on;
PlotPSTH(time_axis,TP_ALL_Eff_Sacc8_Avg,colorsC(1,:));
PlotPSTH(time_axis,TA_ALL_Eff_Sacc8_Avg,colorsC(2,:));
% Labeling
xlabel('Time (s)');
ylabel('Firing Rate (spikes/s)');
title('Efficient search wirh SaccNum=8');
legend('TP', 'TA');
grid on;
hold off;


%% Plotting for InEff

time_axis = linspace(-0.6, 1.0, 1600);
colorsC = get_distinguishable_colors(2); % You can also use other colormap functions like parula, jet, etc.
 

figure; % Creates a new figure window
 % Create subplot in a 2x4 grid at the i-th position
 %sacc1
subplot(2, 4, 1);
hold on;
PlotPSTH(time_axis,TP_ALL_Ineff_Sacc1_Avg,colorsC(1,:));
PlotPSTH(time_axis,TA_ALL_Ineff_Sacc1_Avg,colorsC(2,:));
% Labeling
xlabel('Time (s)');
ylabel('Firing Rate (spikes/s)');
title('Efficient search wirh SaccNum=1');
legend('TP', 'TA');
grid on;
hold off;

%sacc2
subplot(2, 4, 2);
hold on;
PlotPSTH(time_axis,TP_ALL_Ineff_Sacc2_Avg,colorsC(1,:));
PlotPSTH(time_axis,TA_ALL_Ineff_Sacc2_Avg,colorsC(2,:));
% Labeling
xlabel('Time (s)');
ylabel('Firing Rate (spikes/s)');
title('Efficient search wirh SaccNum=2');
legend('TP', 'TA');
grid on;
hold off;

%sacc3
subplot(2, 4, 3);
hold on;
PlotPSTH(time_axis,TP_ALL_Ineff_Sacc3_Avg,colorsC(1,:));
PlotPSTH(time_axis,TA_ALL_Ineff_Sacc3_Avg,colorsC(2,:));
% Labeling
xlabel('Time (s)');
ylabel('Firing Rate (spikes/s)');
title('Efficient search wirh SaccNum=3');
legend('TP', 'TA');
grid on;
hold off;

%sacc4
subplot(2, 4, 4);
hold on;
PlotPSTH(time_axis,TP_ALL_Ineff_Sacc4_Avg,colorsC(1,:));
PlotPSTH(time_axis,TA_ALL_Ineff_Sacc4_Avg,colorsC(2,:));
% Labeling
xlabel('Time (s)');
ylabel('Firing Rate (spikes/s)');
title('Efficient search wirh SaccNum=4');
legend('TP', 'TA');
grid on;
hold off;

%sacc5
subplot(2, 4, 5);
hold on;
PlotPSTH(time_axis,TP_ALL_Ineff_Sacc5_Avg,colorsC(1,:));
PlotPSTH(time_axis,TA_ALL_Ineff_Sacc5_Avg,colorsC(2,:));
% Labeling
xlabel('Time (s)');
ylabel('Firing Rate (spikes/s)');
title('Efficient search wirh SaccNum=5');
legend('TP', 'TA');
grid on;
hold off;

%sacc6
subplot(2, 4, 6);
hold on;
PlotPSTH(time_axis,TP_ALL_Ineff_Sacc6_Avg,colorsC(1,:));
PlotPSTH(time_axis,TA_ALL_Ineff_Sacc6_Avg,colorsC(2,:));
% Labeling
xlabel('Time (s)');
ylabel('Firing Rate (spikes/s)');
title('Efficient search wirh SaccNum=6');
legend('TP', 'TA');
grid on;
hold off;

%sacc7
subplot(2, 4, 7);
hold on;
PlotPSTH(time_axis,TP_ALL_Ineff_Sacc7_Avg,colorsC(1,:));
PlotPSTH(time_axis,TA_ALL_Ineff_Sacc7_Avg,colorsC(2,:));
% Labeling
xlabel('Time (s)');
ylabel('Firing Rate (spikes/s)');
title('Efficient search wirh SaccNum=7');
legend('TP', 'TA');
grid on;
hold off;

%sacc8
subplot(2, 4, 8);
hold on;
PlotPSTH(time_axis,TP_ALL_Ineff_Sacc8_Avg,colorsC(1,:));
PlotPSTH(time_axis,TA_ALL_Ineff_Sacc8_Avg,colorsC(2,:));
% Labeling
xlabel('Time (s)');
ylabel('Firing Rate (spikes/s)');
title('Efficient search wirh SaccNum=8');
legend('TP', 'TA');
grid on;
hold off;

%% Correct Incorrect separation
% Corr/Incorr Separating
Correct=find(table.Accuracy==1);
InCorrect=find(table.Accuracy==0);

% 1 Saccade Correct/Incorrect
TP_ALL_Eff_Sacc1_temp=intersect(intersect(TempTP,Eff),Sacc1);
TP_ALL_Eff_Sacc1_temp_correct=intersect(TP_ALL_Eff_Sacc1_temp,Correct);
TP_ALL_Eff_Sacc1_temp_Incorrect=intersect(TP_ALL_Eff_Sacc1_temp,InCorrect);
TP_ALL_Eff_Sacc1_Incorrect=table(TP_ALL_Eff_Sacc1_temp_Incorrect,:);
TP_ALL_Eff_Sacc1_correct=table(TP_ALL_Eff_Sacc1_temp_correct,:);
TP_ALL_Eff_Sacc1_Incorrect_Bins=TP_ALL_Eff_Sacc1_Incorrect{:, column_names(1:1600)};
TP_ALL_Eff_Sacc1_correct_Bins=TP_ALL_Eff_Sacc1_correct{:, column_names(1:1600)};
TP_ALL_Eff_Sacc1_Incorrect_Bins_Avg=nanmean(TP_ALL_Eff_Sacc1_Incorrect_Bins,1);%
TP_ALL_Eff_Sacc1_correct_Bins_Avg=nanmean(TP_ALL_Eff_Sacc1_correct_Bins,1);%

TP_ALL_Ineff_Sacc1_temp=intersect(intersect(TempTP,Ineff),Sacc1);
TP_ALL_Ineff_Sacc1_temp_correct=intersect(TP_ALL_Ineff_Sacc1_temp,Correct);
TP_ALL_Ineff_Sacc1_temp_Incorrect=intersect(TP_ALL_Ineff_Sacc1_temp,InCorrect);
TP_ALL_Ineff_Sacc1_Incorrect=table(TP_ALL_Ineff_Sacc1_temp_Incorrect,:);
TP_ALL_Ineff_Sacc1_correct=table(TP_ALL_Ineff_Sacc1_temp_correct,:);
TP_ALL_Ineff_Sacc1_Incorrect_Bins=TP_ALL_Ineff_Sacc1_Incorrect{:, column_names(1:1600)};
TP_ALL_Ineff_Sacc1_correct_Bins=TP_ALL_Ineff_Sacc1_correct{:, column_names(1:1600)};
TP_ALL_Ineff_Sacc1_Incorrect_Bins_Avg=nanmean(TP_ALL_Ineff_Sacc1_Incorrect_Bins,1);%
TP_ALL_Ineff_Sacc1_correct_Bins_Avg=nanmean(TP_ALL_Ineff_Sacc1_correct_Bins,1);%

TA_ALL_Eff_Sacc1_temp=intersect(intersect(TempTA,Eff),Sacc1);
TA_ALL_Eff_Sacc1_temp_correct=intersect(TA_ALL_Eff_Sacc1_temp,Correct);
TA_ALL_Eff_Sacc1_temp_Incorrect=intersect(TA_ALL_Eff_Sacc1_temp,InCorrect);
TA_ALL_Eff_Sacc1_Incorrect=table(TA_ALL_Eff_Sacc1_temp_Incorrect,:);
TA_ALL_Eff_Sacc1_correct=table(TA_ALL_Eff_Sacc1_temp_correct,:);
TA_ALL_Eff_Sacc1_Incorrect_Bins=TA_ALL_Eff_Sacc1_Incorrect{:, column_names(1:1600)};
TA_ALL_Eff_Sacc1_correct_Bins=TA_ALL_Eff_Sacc1_correct{:, column_names(1:1600)};
TA_ALL_Eff_Sacc1_Incorrect_Bins_Avg=nanmean(TA_ALL_Eff_Sacc1_Incorrect_Bins,1);%
TA_ALL_Eff_Sacc1_correct_Bins_Avg=nanmean(TA_ALL_Eff_Sacc1_correct_Bins,1);%

TA_ALL_Ineff_Sacc1_temp=intersect(intersect(TempTA,Ineff),Sacc1);
TA_ALL_Ineff_Sacc1_temp_correct=intersect(TA_ALL_Ineff_Sacc1_temp,Correct);
TA_ALL_Ineff_Sacc1_temp_Incorrect=intersect(TA_ALL_Ineff_Sacc1_temp,InCorrect);
TA_ALL_Ineff_Sacc1_Incorrect=table(TA_ALL_Ineff_Sacc1_temp_Incorrect,:);
TA_ALL_Ineff_Sacc1_correct=table(TA_ALL_Ineff_Sacc1_temp_correct,:);
TA_ALL_Ineff_Sacc1_Incorrect_Bins=TA_ALL_Ineff_Sacc1_Incorrect{:, column_names(1:1600)};
TA_ALL_Ineff_Sacc1_correct_Bins=TA_ALL_Ineff_Sacc1_correct{:, column_names(1:1600)};
TA_ALL_Ineff_Sacc1_Incorrect_Bins_Avg=nanmean(TA_ALL_Ineff_Sacc1_Incorrect_Bins,1);%
TA_ALL_Ineff_Sacc1_correct_Bins_Avg=nanmean(TA_ALL_Ineff_Sacc1_correct_Bins,1);%



% 2 Saccade Correct/Incorrect
TP_ALL_Eff_Sacc2_temp=intersect(intersect(TempTP,Eff),Sacc2);
TP_ALL_Eff_Sacc2_temp_correct=intersect(TP_ALL_Eff_Sacc2_temp,Correct);
TP_ALL_Eff_Sacc2_temp_Incorrect=intersect(TP_ALL_Eff_Sacc2_temp,InCorrect);
TP_ALL_Eff_Sacc2_Incorrect=table(TP_ALL_Eff_Sacc2_temp_Incorrect,:);
TP_ALL_Eff_Sacc2_correct=table(TP_ALL_Eff_Sacc2_temp_correct,:);
TP_ALL_Eff_Sacc2_Incorrect_Bins=TP_ALL_Eff_Sacc2_Incorrect{:, column_names(1:1600)};
TP_ALL_Eff_Sacc2_correct_Bins=TP_ALL_Eff_Sacc2_correct{:, column_names(1:1600)};
TP_ALL_Eff_Sacc2_Incorrect_Bins_Avg=nanmean(TP_ALL_Eff_Sacc2_Incorrect_Bins,1);%
TP_ALL_Eff_Sacc2_correct_Bins_Avg=nanmean(TP_ALL_Eff_Sacc2_correct_Bins,1);%

TP_ALL_Ineff_Sacc2_temp=intersect(intersect(TempTP,Ineff),Sacc2);
TP_ALL_Ineff_Sacc2_temp_correct=intersect(TP_ALL_Ineff_Sacc2_temp,Correct);
TP_ALL_Ineff_Sacc2_temp_Incorrect=intersect(TP_ALL_Ineff_Sacc2_temp,InCorrect);
TP_ALL_Ineff_Sacc2_Incorrect=table(TP_ALL_Ineff_Sacc2_temp_Incorrect,:);
TP_ALL_Ineff_Sacc2_correct=table(TP_ALL_Ineff_Sacc2_temp_correct,:);
TP_ALL_Ineff_Sacc2_Incorrect_Bins=TP_ALL_Ineff_Sacc2_Incorrect{:, column_names(1:1600)};
TP_ALL_Ineff_Sacc2_correct_Bins=TP_ALL_Ineff_Sacc2_correct{:, column_names(1:1600)};
TP_ALL_Ineff_Sacc2_Incorrect_Bins_Avg=nanmean(TP_ALL_Ineff_Sacc2_Incorrect_Bins,1);%
TP_ALL_Ineff_Sacc2_correct_Bins_Avg=nanmean(TP_ALL_Ineff_Sacc2_correct_Bins,1);%

TA_ALL_Eff_Sacc2_temp=intersect(intersect(TempTA,Eff),Sacc2);
TA_ALL_Eff_Sacc2_temp_correct=intersect(TA_ALL_Eff_Sacc2_temp,Correct);
TA_ALL_Eff_Sacc2_temp_Incorrect=intersect(TA_ALL_Eff_Sacc2_temp,InCorrect);
TA_ALL_Eff_Sacc2_Incorrect=table(TA_ALL_Eff_Sacc2_temp_Incorrect,:);
TA_ALL_Eff_Sacc2_correct=table(TA_ALL_Eff_Sacc2_temp_correct,:);
TA_ALL_Eff_Sacc2_Incorrect_Bins=TA_ALL_Eff_Sacc2_Incorrect{:, column_names(1:1600)};
TA_ALL_Eff_Sacc2_correct_Bins=TA_ALL_Eff_Sacc2_correct{:, column_names(1:1600)};
TA_ALL_Eff_Sacc2_Incorrect_Bins_Avg=nanmean(TA_ALL_Eff_Sacc2_Incorrect_Bins,1);%
TA_ALL_Eff_Sacc2_correct_Bins_Avg=nanmean(TA_ALL_Eff_Sacc2_correct_Bins,1);%

TA_ALL_Ineff_Sacc2_temp=intersect(intersect(TempTA,Ineff),Sacc2);
TA_ALL_Ineff_Sacc2_temp_correct=intersect(TA_ALL_Ineff_Sacc2_temp,Correct);
TA_ALL_Ineff_Sacc2_temp_Incorrect=intersect(TA_ALL_Ineff_Sacc2_temp,InCorrect);
TA_ALL_Ineff_Sacc2_Incorrect=table(TA_ALL_Ineff_Sacc2_temp_Incorrect,:);
TA_ALL_Ineff_Sacc2_correct=table(TA_ALL_Ineff_Sacc2_temp_correct,:);
TA_ALL_Ineff_Sacc2_Incorrect_Bins=TA_ALL_Ineff_Sacc2_Incorrect{:, column_names(1:1600)};
TA_ALL_Ineff_Sacc2_correct_Bins=TA_ALL_Ineff_Sacc2_correct{:, column_names(1:1600)};
TA_ALL_Ineff_Sacc2_Incorrect_Bins_Avg=nanmean(TA_ALL_Ineff_Sacc2_Incorrect_Bins,1);%
TA_ALL_Ineff_Sacc2_correct_Bins_Avg=nanmean(TA_ALL_Ineff_Sacc2_correct_Bins,1);%

% 3 Saccade Correct/Incorrect
TP_ALL_Eff_Sacc3_temp=intersect(intersect(TempTP,Eff),Sacc3);
TP_ALL_Eff_Sacc3_temp_correct=intersect(TP_ALL_Eff_Sacc3_temp,Correct);
TP_ALL_Eff_Sacc3_temp_Incorrect=intersect(TP_ALL_Eff_Sacc3_temp,InCorrect);
TP_ALL_Eff_Sacc3_Incorrect=table(TP_ALL_Eff_Sacc3_temp_Incorrect,:);
TP_ALL_Eff_Sacc3_correct=table(TP_ALL_Eff_Sacc3_temp_correct,:);
TP_ALL_Eff_Sacc3_Incorrect_Bins=TP_ALL_Eff_Sacc3_Incorrect{:, column_names(1:1600)};
TP_ALL_Eff_Sacc3_correct_Bins=TP_ALL_Eff_Sacc3_correct{:, column_names(1:1600)};
TP_ALL_Eff_Sacc3_Incorrect_Bins_Avg=nanmean(TP_ALL_Eff_Sacc3_Incorrect_Bins,1);%
TP_ALL_Eff_Sacc3_correct_Bins_Avg=nanmean(TP_ALL_Eff_Sacc3_correct_Bins,1);%

TP_ALL_Ineff_Sacc3_temp=intersect(intersect(TempTP,Ineff),Sacc3);
TP_ALL_Ineff_Sacc3_temp_correct=intersect(TP_ALL_Ineff_Sacc3_temp,Correct);
TP_ALL_Ineff_Sacc3_temp_Incorrect=intersect(TP_ALL_Ineff_Sacc3_temp,InCorrect);
TP_ALL_Ineff_Sacc3_Incorrect=table(TP_ALL_Ineff_Sacc3_temp_Incorrect,:);
TP_ALL_Ineff_Sacc3_correct=table(TP_ALL_Ineff_Sacc3_temp_correct,:);
TP_ALL_Ineff_Sacc3_Incorrect_Bins=TP_ALL_Ineff_Sacc3_Incorrect{:, column_names(1:1600)};
TP_ALL_Ineff_Sacc3_correct_Bins=TP_ALL_Ineff_Sacc3_correct{:, column_names(1:1600)};
TP_ALL_Ineff_Sacc3_Incorrect_Bins_Avg=nanmean(TP_ALL_Ineff_Sacc3_Incorrect_Bins,1);%
TP_ALL_Ineff_Sacc3_correct_Bins_Avg=nanmean(TP_ALL_Ineff_Sacc3_correct_Bins,1);%

TA_ALL_Eff_Sacc3_temp=intersect(intersect(TempTA,Eff),Sacc3);
TA_ALL_Eff_Sacc3_temp_correct=intersect(TA_ALL_Eff_Sacc3_temp,Correct);
TA_ALL_Eff_Sacc3_temp_Incorrect=intersect(TA_ALL_Eff_Sacc3_temp,InCorrect);
TA_ALL_Eff_Sacc3_Incorrect=table(TA_ALL_Eff_Sacc3_temp_Incorrect,:);
TA_ALL_Eff_Sacc3_correct=table(TA_ALL_Eff_Sacc3_temp_correct,:);
TA_ALL_Eff_Sacc3_Incorrect_Bins=TA_ALL_Eff_Sacc3_Incorrect{:, column_names(1:1600)};
TA_ALL_Eff_Sacc3_correct_Bins=TA_ALL_Eff_Sacc3_correct{:, column_names(1:1600)};
TA_ALL_Eff_Sacc3_Incorrect_Bins_Avg=nanmean(TA_ALL_Eff_Sacc3_Incorrect_Bins,1);%
TA_ALL_Eff_Sacc3_correct_Bins_Avg=nanmean(TA_ALL_Eff_Sacc3_correct_Bins,1);%

TA_ALL_Ineff_Sacc3_temp=intersect(intersect(TempTA,Ineff),Sacc3);
TA_ALL_Ineff_Sacc3_temp_correct=intersect(TA_ALL_Ineff_Sacc3_temp,Correct);
TA_ALL_Ineff_Sacc3_temp_Incorrect=intersect(TA_ALL_Ineff_Sacc3_temp,InCorrect);
TA_ALL_Ineff_Sacc3_Incorrect=table(TA_ALL_Ineff_Sacc3_temp_Incorrect,:);
TA_ALL_Ineff_Sacc3_correct=table(TA_ALL_Ineff_Sacc3_temp_correct,:);
TA_ALL_Ineff_Sacc3_Incorrect_Bins=TA_ALL_Ineff_Sacc3_Incorrect{:, column_names(1:1600)};
TA_ALL_Ineff_Sacc3_correct_Bins=TA_ALL_Ineff_Sacc3_correct{:, column_names(1:1600)};
TA_ALL_Ineff_Sacc3_Incorrect_Bins_Avg=nanmean(TA_ALL_Ineff_Sacc3_Incorrect_Bins,1);%
TA_ALL_Ineff_Sacc3_correct_Bins_Avg=nanmean(TA_ALL_Ineff_Sacc3_correct_Bins,1);%

% 4 Saccade Correct/Incorrect
TP_ALL_Eff_Sacc4_temp=intersect(intersect(TempTP,Eff),Sacc4);
TP_ALL_Eff_Sacc4_temp_correct=intersect(TP_ALL_Eff_Sacc4_temp,Correct);
TP_ALL_Eff_Sacc4_temp_Incorrect=intersect(TP_ALL_Eff_Sacc4_temp,InCorrect);
TP_ALL_Eff_Sacc4_Incorrect=table(TP_ALL_Eff_Sacc4_temp_Incorrect,:);
TP_ALL_Eff_Sacc4_correct=table(TP_ALL_Eff_Sacc4_temp_correct,:);
TP_ALL_Eff_Sacc4_Incorrect_Bins=TP_ALL_Eff_Sacc4_Incorrect{:, column_names(1:1600)};
TP_ALL_Eff_Sacc4_correct_Bins=TP_ALL_Eff_Sacc4_correct{:, column_names(1:1600)};
TP_ALL_Eff_Sacc4_Incorrect_Bins_Avg=nanmean(TP_ALL_Eff_Sacc4_Incorrect_Bins,1);%
TP_ALL_Eff_Sacc4_correct_Bins_Avg=nanmean(TP_ALL_Eff_Sacc4_correct_Bins,1);%

TP_ALL_Ineff_Sacc4_temp=intersect(intersect(TempTP,Ineff),Sacc4);
TP_ALL_Ineff_Sacc4_temp_correct=intersect(TP_ALL_Ineff_Sacc4_temp,Correct);
TP_ALL_Ineff_Sacc4_temp_Incorrect=intersect(TP_ALL_Ineff_Sacc4_temp,InCorrect);
TP_ALL_Ineff_Sacc4_Incorrect=table(TP_ALL_Ineff_Sacc4_temp_Incorrect,:);
TP_ALL_Ineff_Sacc4_correct=table(TP_ALL_Ineff_Sacc4_temp_correct,:);
TP_ALL_Ineff_Sacc4_Incorrect_Bins=TP_ALL_Ineff_Sacc4_Incorrect{:, column_names(1:1600)};
TP_ALL_Ineff_Sacc4_correct_Bins=TP_ALL_Ineff_Sacc4_correct{:, column_names(1:1600)};
TP_ALL_Ineff_Sacc4_Incorrect_Bins_Avg=nanmean(TP_ALL_Ineff_Sacc4_Incorrect_Bins,1);%
TP_ALL_Ineff_Sacc4_correct_Bins_Avg=nanmean(TP_ALL_Ineff_Sacc4_correct_Bins,1);%

TA_ALL_Eff_Sacc4_temp=intersect(intersect(TempTA,Eff),Sacc4);
TA_ALL_Eff_Sacc4_temp_correct=intersect(TA_ALL_Eff_Sacc4_temp,Correct);
TA_ALL_Eff_Sacc4_temp_Incorrect=intersect(TA_ALL_Eff_Sacc4_temp,InCorrect);
TA_ALL_Eff_Sacc4_Incorrect=table(TA_ALL_Eff_Sacc4_temp_Incorrect,:);
TA_ALL_Eff_Sacc4_correct=table(TA_ALL_Eff_Sacc4_temp_correct,:);
TA_ALL_Eff_Sacc4_Incorrect_Bins=TA_ALL_Eff_Sacc4_Incorrect{:, column_names(1:1600)};
TA_ALL_Eff_Sacc4_correct_Bins=TA_ALL_Eff_Sacc4_correct{:, column_names(1:1600)};
TA_ALL_Eff_Sacc4_Incorrect_Bins_Avg=nanmean(TA_ALL_Eff_Sacc4_Incorrect_Bins,1);%
TA_ALL_Eff_Sacc4_correct_Bins_Avg=nanmean(TA_ALL_Eff_Sacc4_correct_Bins,1);%

TA_ALL_Ineff_Sacc4_temp=intersect(intersect(TempTA,Ineff),Sacc4);
TA_ALL_Ineff_Sacc4_temp_correct=intersect(TA_ALL_Ineff_Sacc4_temp,Correct);
TA_ALL_Ineff_Sacc4_temp_Incorrect=intersect(TA_ALL_Ineff_Sacc4_temp,InCorrect);
TA_ALL_Ineff_Sacc4_Incorrect=table(TA_ALL_Ineff_Sacc4_temp_Incorrect,:);
TA_ALL_Ineff_Sacc4_correct=table(TA_ALL_Ineff_Sacc4_temp_correct,:);
TA_ALL_Ineff_Sacc4_Incorrect_Bins=TA_ALL_Ineff_Sacc4_Incorrect{:, column_names(1:1600)};
TA_ALL_Ineff_Sacc4_correct_Bins=TA_ALL_Ineff_Sacc4_correct{:, column_names(1:1600)};
TA_ALL_Ineff_Sacc4_Incorrect_Bins_Avg=nanmean(TA_ALL_Ineff_Sacc4_Incorrect_Bins,1);%
TA_ALL_Ineff_Sacc4_correct_Bins_Avg=nanmean(TA_ALL_Ineff_Sacc4_correct_Bins,1);%
%% Plotting for Eff Correct/Incorrect

time_axis = linspace(-0.6, 1.0, 1600);
colorsC = get_distinguishable_colors(2); % You can also use other colormap functions like parula, jet, etc.
 

figure; % Creates a new figure window
 % Create subplot in a 2x4 grid at the i-th position
 %sacc1
subplot(2, 4, 1);
hold on;
PlotPSTH(time_axis,TP_ALL_Eff_Sacc1_correct_Bins_Avg,colorsC(1,:));
PlotPSTH(time_axis,TA_ALL_Eff_Sacc1_correct_Bins_Avg,colorsC(2,:));
% Labeling
xlabel('Time (s)');
ylabel('Firing Rate (spikes/s)');
title('Efficient search with SaccNum=1 and Correct');
legend('TP', 'TA');
grid on;
hold off;


%sacc2
subplot(2, 4, 2);
hold on;
PlotPSTH(time_axis,TP_ALL_Eff_Sacc2_correct_Bins_Avg,colorsC(1,:));
PlotPSTH(time_axis,TA_ALL_Eff_Sacc2_correct_Bins_Avg,colorsC(2,:));
% Labeling
xlabel('Time (s)');
ylabel('Firing Rate (spikes/s)');
title('Efficient search with SaccNum=2 and Correct');
legend('TP', 'TA');
grid on;
hold off;


%sacc3
subplot(2, 4, 3);
hold on;
PlotPSTH(time_axis,TP_ALL_Eff_Sacc3_correct_Bins_Avg,colorsC(1,:));
PlotPSTH(time_axis,TA_ALL_Eff_Sacc3_correct_Bins_Avg,colorsC(2,:));
% Labeling
xlabel('Time (s)');
ylabel('Firing Rate (spikes/s)');
title('Efficient search with SaccNum=3 and Correct');
legend('TP', 'TA');
grid on;
hold off;

%sacc4
subplot(2, 4, 4);
hold on;
PlotPSTH(time_axis,TP_ALL_Eff_Sacc4_correct_Bins_Avg,colorsC(1,:));
PlotPSTH(time_axis,TA_ALL_Eff_Sacc4_correct_Bins_Avg,colorsC(2,:));
% Labeling
xlabel('Time (s)');
ylabel('Firing Rate (spikes/s)');
title('Efficient search with SaccNum=4 and Correct');
legend('TP', 'TA');
grid on;
hold off;

%sacc1 incorr
subplot(2, 4, 5);
hold on;
PlotPSTH(time_axis,TP_ALL_Eff_Sacc1_Incorrect_Bins_Avg,colorsC(1,:));
PlotPSTH(time_axis,TA_ALL_Eff_Sacc1_Incorrect_Bins_Avg,colorsC(2,:));
% Labeling
xlabel('Time (s)');
ylabel('Firing Rate (spikes/s)');
title('Efficient search with SaccNum=1 and Incorrect');
legend('TP', 'TA');
grid on;
hold off;

%sacc2 incorr
subplot(2, 4, 6);
hold on;
PlotPSTH(time_axis,TP_ALL_Eff_Sacc2_Incorrect_Bins_Avg,colorsC(1,:));
PlotPSTH(time_axis,TA_ALL_Eff_Sacc2_Incorrect_Bins_Avg,colorsC(2,:));
% Labeling
xlabel('Time (s)');
ylabel('Firing Rate (spikes/s)');
title('Efficient search with SaccNum=2 and Incorrect');
legend('TP', 'TA');
grid on;
hold off;

%sacc3 incorr
subplot(2, 4, 7);
hold on;
PlotPSTH(time_axis,TP_ALL_Eff_Sacc3_Incorrect_Bins_Avg,colorsC(1,:));
PlotPSTH(time_axis,TA_ALL_Eff_Sacc3_Incorrect_Bins_Avg,colorsC(2,:));
% Labeling
xlabel('Time (s)');
ylabel('Firing Rate (spikes/s)');
title('Efficient search with SaccNum=3 and Incorrect');
legend('TP', 'TA');
grid on;
hold off;

%sacc4 incorr
subplot(2, 4, 8);
hold on;
PlotPSTH(time_axis,TP_ALL_Eff_Sacc4_Incorrect_Bins_Avg,colorsC(1,:));
PlotPSTH(time_axis,TA_ALL_Eff_Sacc4_Incorrect_Bins_Avg,colorsC(2,:));
% Labeling
xlabel('Time (s)');
ylabel('Firing Rate (spikes/s)');
title('Efficient search with SaccNum=4 and Incorrect');
legend('TP', 'TA');
grid on;
hold off;


%% Plotting for InEff Correct/Incorrect

time_axis = linspace(-0.6, 1.0, 1600);
colorsC = get_distinguishable_colors(2); % You can also use other colormap functions like parula, jet, etc.
 

figure; % Creates a new figure window
 % Create subplot in a 2x4 grid at the i-th position
 %sacc1
subplot(2, 4, 1);
hold on;
PlotPSTH(time_axis,TP_ALL_Ineff_Sacc1_correct_Bins_Avg,colorsC(1,:));
PlotPSTH(time_axis,TA_ALL_Ineff_Sacc1_correct_Bins_Avg,colorsC(2,:));
% Labeling
xlabel('Time (s)');
ylabel('Firing Rate (spikes/s)');
title('InEfficient search with SaccNum=1 and Correct');
legend('TP', 'TA');
grid on;
hold off;


%sacc2
subplot(2, 4, 2);
hold on;
PlotPSTH(time_axis,TP_ALL_Ineff_Sacc2_correct_Bins_Avg,colorsC(1,:));
PlotPSTH(time_axis,TA_ALL_Ineff_Sacc2_correct_Bins_Avg,colorsC(2,:));
% Labeling
xlabel('Time (s)');
ylabel('Firing Rate (spikes/s)');
title('InEfficient search with SaccNum=2 and Correct');
legend('TP', 'TA');
grid on;
hold off;


%sacc3
subplot(2, 4, 3);
hold on;
PlotPSTH(time_axis,TP_ALL_Ineff_Sacc3_correct_Bins_Avg,colorsC(1,:));
PlotPSTH(time_axis,TA_ALL_Ineff_Sacc3_correct_Bins_Avg,colorsC(2,:));
% Labeling
xlabel('Time (s)');
ylabel('Firing Rate (spikes/s)');
title('InEfficient search with SaccNum=3 and Correct');
legend('TP', 'TA');
grid on;
hold off;

%sacc4
subplot(2, 4, 4);
hold on;
PlotPSTH(time_axis,TP_ALL_Ineff_Sacc4_correct_Bins_Avg,colorsC(1,:));
PlotPSTH(time_axis,TA_ALL_Ineff_Sacc4_correct_Bins_Avg,colorsC(2,:));
% Labeling
xlabel('Time (s)');
ylabel('Firing Rate (spikes/s)');
title('InEfficient search with SaccNum=4 and Correct');
legend('TP', 'TA');
grid on;
hold off;

%sacc1 incorr
subplot(2, 4, 5);
hold on;
PlotPSTH(time_axis,TP_ALL_Ineff_Sacc1_Incorrect_Bins_Avg,colorsC(1,:));
PlotPSTH(time_axis,TA_ALL_Ineff_Sacc1_Incorrect_Bins_Avg,colorsC(2,:));
% Labeling
xlabel('Time (s)');
ylabel('Firing Rate (spikes/s)');
title('InEfficient search with SaccNum=1 and Incorrect');
legend('TP', 'TA');
grid on;
hold off;

%sacc2 incorr
subplot(2, 4, 6);
hold on;
PlotPSTH(time_axis,TP_ALL_Ineff_Sacc2_Incorrect_Bins_Avg,colorsC(1,:));
PlotPSTH(time_axis,TA_ALL_Ineff_Sacc2_Incorrect_Bins_Avg,colorsC(2,:));
% Labeling
xlabel('Time (s)');
ylabel('Firing Rate (spikes/s)');
title('InEfficient search with SaccNum=2 and Incorrect');
legend('TP', 'TA');
grid on;
hold off;

%sacc3 incorr
subplot(2, 4, 7);
hold on;
PlotPSTH(time_axis,TP_ALL_Ineff_Sacc3_Incorrect_Bins_Avg,colorsC(1,:));
PlotPSTH(time_axis,TA_ALL_Ineff_Sacc3_Incorrect_Bins_Avg,colorsC(2,:));
% Labeling
xlabel('Time (s)');
ylabel('Firing Rate (spikes/s)');
title('InEfficient search with SaccNum=3 and Incorrect');
legend('TP', 'TA');
grid on;
hold off;

%sacc4 incorr
subplot(2, 4, 8);
hold on;
PlotPSTH(time_axis,TP_ALL_Ineff_Sacc4_Incorrect_Bins_Avg,colorsC(1,:));
PlotPSTH(time_axis,TA_ALL_Ineff_Sacc4_Incorrect_Bins_Avg,colorsC(2,:));
% Labeling
xlabel('Time (s)');
ylabel('Firing Rate (spikes/s)');
title('InEfficient search with SaccNum=4 and Incorrect');
legend('TP', 'TA');
grid on;
hold off;

%% Separating NoSacc(Reject) / 1 Sacc / 2 Sacc / >2 Sacc







%% Functions
% Plot PSTH
% Gives mean bins of the data 
function PlotPSTH(time_axis,meanbins,color)

sigma1 = 4; % Adjust sigma value as needed
 % Adjust window size as needed
smoothed_data = imgaussfilt(meanbins, sigma1);

% Plot PSTH as a continuous line

plot(time_axis, smoothed_data, 'Color', color, 'LineWidth', 2);

end


function colors = get_distinguishable_colors(n)
    % Generate a set of distinguishable colors
    golden_ratio_conjugate = (1 + sqrt(5)) / 2;
    hue = mod((0:n-1)' / golden_ratio_conjugate, 1);
    saturation = 0.6;
    value = 0.9;
    colors = hsv2rgb([hue, saturation * ones(n, 1), value * ones(n, 1)]);
end






