% Sacc analysis for SNr

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

% Corr/Incorr Separating
Correct=find(table.Accuracy==1);
InCorrect=find(table.Accuracy==0);

% Sacc Numbers separation
Reject=find(isnan(table.SaccNum));
Sacc1=find(table.SaccNum==1);
Sacc2=find(table.SaccNum==2);
SaccMore2=find(table.SaccNum>=2);

% Reject Correct/Incorrect
TP_ALL_Eff_Reject_temp=intersect(intersect(TempTP,Eff),Reject);
TP_ALL_Eff_Reject_temp_correct=intersect(TP_ALL_Eff_Reject_temp,Correct);
TP_ALL_Eff_Reject_temp_Incorrect=intersect(TP_ALL_Eff_Reject_temp,InCorrect);
TP_ALL_Eff_Reject_Incorrect=table(TP_ALL_Eff_Reject_temp_Incorrect,:);
TP_ALL_Eff_Reject_correct=table(TP_ALL_Eff_Reject_temp_correct,:);
TP_ALL_Eff_Reject_Incorrect_Bins=TP_ALL_Eff_Reject_Incorrect{:, column_names(1:1600)};
TP_ALL_Eff_Reject_correct_Bins=TP_ALL_Eff_Reject_correct{:, column_names(1:1600)};
TP_ALL_Eff_Reject_Incorrect_Bins_Avg=nanmean(TP_ALL_Eff_Reject_Incorrect_Bins,1);%
TP_ALL_Eff_Reject_correct_Bins_Avg=nanmean(TP_ALL_Eff_Reject_correct_Bins,1);%

TP_ALL_Ineff_Reject_temp=intersect(intersect(TempTP,Ineff),Reject);
TP_ALL_Ineff_Reject_temp_correct=intersect(TP_ALL_Ineff_Reject_temp,Correct);
TP_ALL_Ineff_Reject_temp_Incorrect=intersect(TP_ALL_Ineff_Reject_temp,InCorrect);
TP_ALL_Ineff_Reject_Incorrect=table(TP_ALL_Ineff_Reject_temp_Incorrect,:);
TP_ALL_Ineff_Reject_correct=table(TP_ALL_Ineff_Reject_temp_correct,:);
TP_ALL_Ineff_Reject_Incorrect_Bins=TP_ALL_Ineff_Reject_Incorrect{:, column_names(1:1600)};
TP_ALL_Ineff_Reject_correct_Bins=TP_ALL_Ineff_Reject_correct{:, column_names(1:1600)};
TP_ALL_Ineff_Reject_Incorrect_Bins_Avg=nanmean(TP_ALL_Ineff_Reject_Incorrect_Bins,1);%
TP_ALL_Ineff_Reject_correct_Bins_Avg=nanmean(TP_ALL_Ineff_Reject_correct_Bins,1);%

TA_ALL_Eff_Reject_temp=intersect(intersect(TempTA,Eff),Reject);
TA_ALL_Eff_Reject_temp_correct=intersect(TA_ALL_Eff_Reject_temp,Correct);
TA_ALL_Eff_Reject_temp_Incorrect=intersect(TA_ALL_Eff_Reject_temp,InCorrect);
TA_ALL_Eff_Reject_Incorrect=table(TA_ALL_Eff_Reject_temp_Incorrect,:);
TA_ALL_Eff_Reject_correct=table(TA_ALL_Eff_Reject_temp_correct,:);
TA_ALL_Eff_Reject_Incorrect_Bins=TA_ALL_Eff_Reject_Incorrect{:, column_names(1:1600)};
TA_ALL_Eff_Reject_correct_Bins=TA_ALL_Eff_Reject_correct{:, column_names(1:1600)};
TA_ALL_Eff_Reject_Incorrect_Bins_Avg=nanmean(TA_ALL_Eff_Reject_Incorrect_Bins,1);%
TA_ALL_Eff_Reject_correct_Bins_Avg=nanmean(TA_ALL_Eff_Reject_correct_Bins,1);%

TA_ALL_Ineff_Reject_temp=intersect(intersect(TempTA,Ineff),Reject);
TA_ALL_Ineff_Reject_temp_correct=intersect(TA_ALL_Ineff_Reject_temp,Correct);
TA_ALL_Ineff_Reject_temp_Incorrect=intersect(TA_ALL_Ineff_Reject_temp,InCorrect);
TA_ALL_Ineff_Reject_Incorrect=table(TA_ALL_Ineff_Reject_temp_Incorrect,:);
TA_ALL_Ineff_Reject_correct=table(TA_ALL_Ineff_Reject_temp_correct,:);
TA_ALL_Ineff_Reject_Incorrect_Bins=TA_ALL_Ineff_Reject_Incorrect{:, column_names(1:1600)};
TA_ALL_Ineff_Reject_correct_Bins=TA_ALL_Ineff_Reject_correct{:, column_names(1:1600)};
TA_ALL_Ineff_Reject_Incorrect_Bins_Avg=nanmean(TA_ALL_Ineff_Reject_Incorrect_Bins,1);%
TA_ALL_Ineff_Reject_correct_Bins_Avg=nanmean(TA_ALL_Ineff_Reject_correct_Bins,1);%

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



% More than 2 Saccade Correct/Incorrect
TP_ALL_Eff_SaccMore2_temp=intersect(intersect(TempTP,Eff),SaccMore2);
TP_ALL_Eff_SaccMore2_temp_correct=intersect(TP_ALL_Eff_SaccMore2_temp,Correct);
TP_ALL_Eff_SaccMore2_temp_Incorrect=intersect(TP_ALL_Eff_SaccMore2_temp,InCorrect);
TP_ALL_Eff_SaccMore2_Incorrect=table(TP_ALL_Eff_SaccMore2_temp_Incorrect,:);
TP_ALL_Eff_SaccMore2_correct=table(TP_ALL_Eff_SaccMore2_temp_correct,:);
TP_ALL_Eff_SaccMore2_Incorrect_Bins=TP_ALL_Eff_SaccMore2_Incorrect{:, column_names(1:1600)};
TP_ALL_Eff_SaccMore2_correct_Bins=TP_ALL_Eff_SaccMore2_correct{:, column_names(1:1600)};
TP_ALL_Eff_SaccMore2_Incorrect_Bins_Avg=nanmean(TP_ALL_Eff_SaccMore2_Incorrect_Bins,1);%
TP_ALL_Eff_SaccMore2_correct_Bins_Avg=nanmean(TP_ALL_Eff_SaccMore2_correct_Bins,1);%

TP_ALL_Ineff_SaccMore2_temp=intersect(intersect(TempTP,Ineff),SaccMore2);
TP_ALL_Ineff_SaccMore2_temp_correct=intersect(TP_ALL_Ineff_SaccMore2_temp,Correct);
TP_ALL_Ineff_SaccMore2_temp_Incorrect=intersect(TP_ALL_Ineff_SaccMore2_temp,InCorrect);
TP_ALL_Ineff_SaccMore2_Incorrect=table(TP_ALL_Ineff_SaccMore2_temp_Incorrect,:);
TP_ALL_Ineff_SaccMore2_correct=table(TP_ALL_Ineff_SaccMore2_temp_correct,:);
TP_ALL_Ineff_SaccMore2_Incorrect_Bins=TP_ALL_Ineff_SaccMore2_Incorrect{:, column_names(1:1600)};
TP_ALL_Ineff_SaccMore2_correct_Bins=TP_ALL_Ineff_SaccMore2_correct{:, column_names(1:1600)};
TP_ALL_Ineff_SaccMore2_Incorrect_Bins_Avg=nanmean(TP_ALL_Ineff_SaccMore2_Incorrect_Bins,1);%
TP_ALL_Ineff_SaccMore2_correct_Bins_Avg=nanmean(TP_ALL_Ineff_SaccMore2_correct_Bins,1);%

TA_ALL_Eff_SaccMore2_temp=intersect(intersect(TempTA,Eff),SaccMore2);
TA_ALL_Eff_SaccMore2_temp_correct=intersect(TA_ALL_Eff_SaccMore2_temp,Correct);
TA_ALL_Eff_SaccMore2_temp_Incorrect=intersect(TA_ALL_Eff_SaccMore2_temp,InCorrect);
TA_ALL_Eff_SaccMore2_Incorrect=table(TA_ALL_Eff_SaccMore2_temp_Incorrect,:);
TA_ALL_Eff_SaccMore2_correct=table(TA_ALL_Eff_SaccMore2_temp_correct,:);
TA_ALL_Eff_SaccMore2_Incorrect_Bins=TA_ALL_Eff_SaccMore2_Incorrect{:, column_names(1:1600)};
TA_ALL_Eff_SaccMore2_correct_Bins=TA_ALL_Eff_SaccMore2_correct{:, column_names(1:1600)};
TA_ALL_Eff_SaccMore2_Incorrect_Bins_Avg=nanmean(TA_ALL_Eff_SaccMore2_Incorrect_Bins,1);%
TA_ALL_Eff_SaccMore2_correct_Bins_Avg=nanmean(TA_ALL_Eff_SaccMore2_correct_Bins,1);%

TA_ALL_Ineff_SaccMore2_temp=intersect(intersect(TempTA,Ineff),SaccMore2);
TA_ALL_Ineff_SaccMore2_temp_correct=intersect(TA_ALL_Ineff_SaccMore2_temp,Correct);
TA_ALL_Ineff_SaccMore2_temp_Incorrect=intersect(TA_ALL_Ineff_SaccMore2_temp,InCorrect);
TA_ALL_Ineff_SaccMore2_Incorrect=table(TA_ALL_Ineff_SaccMore2_temp_Incorrect,:);
TA_ALL_Ineff_SaccMore2_correct=table(TA_ALL_Ineff_SaccMore2_temp_correct,:);
TA_ALL_Ineff_SaccMore2_Incorrect_Bins=TA_ALL_Ineff_SaccMore2_Incorrect{:, column_names(1:1600)};
TA_ALL_Ineff_SaccMore2_correct_Bins=TA_ALL_Ineff_SaccMore2_correct{:, column_names(1:1600)};
TA_ALL_Ineff_SaccMore2_Incorrect_Bins_Avg=nanmean(TA_ALL_Ineff_SaccMore2_Incorrect_Bins,1);%
TA_ALL_Ineff_SaccMore2_correct_Bins_Avg=nanmean(TA_ALL_Ineff_SaccMore2_correct_Bins,1);%


%% Plotting for Eff Correct/Incorrect

time_axis = linspace(-0.6, 1.0, 1600);
colorsC = get_distinguishable_colors(2); % You can also use other colormap functions like parula, jet, etc.
 

figure; % Creates a new figure window
 % Create subplot in a 2x4 grid at the i-th position
 %Rej
subplot(2, 4, 1);
hold on;
PlotPSTH(time_axis,TP_ALL_Eff_Reject_correct_Bins_Avg,colorsC(1,:));
PlotPSTH(time_axis,TA_ALL_Eff_Reject_correct_Bins_Avg,colorsC(2,:));
% Labeling
xlabel('Time (s)');
ylabel('Firing Rate (spikes/s)');
title('Efficient search with Reject and Correct');
legend('TP', 'TA');
grid on;
hold off;


%sacc1
subplot(2, 4, 2);
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
subplot(2, 4, 3);
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

%saccMore
subplot(2, 4, 4);
hold on;
PlotPSTH(time_axis,TP_ALL_Eff_SaccMore2_correct_Bins_Avg,colorsC(1,:));
PlotPSTH(time_axis,TA_ALL_Eff_SaccMore2_correct_Bins_Avg,colorsC(2,:));
% Labeling
xlabel('Time (s)');
ylabel('Firing Rate (spikes/s)');
title('Efficient search with SaccNum>2 and Correct');
legend('TP', 'TA');
grid on;
hold off;

%Rej incorr
subplot(2, 4, 5);
hold on;
PlotPSTH(time_axis,TP_ALL_Eff_Reject_Incorrect_Bins_Avg,colorsC(1,:));
PlotPSTH(time_axis,TA_ALL_Eff_Reject_Incorrect_Bins_Avg,colorsC(2,:));
% Labeling
xlabel('Time (s)');
ylabel('Firing Rate (spikes/s)');
title('Efficient search with Reject and Incorrect');
legend('TP', 'TA');
grid on;
hold off;

%sacc1 incorr
subplot(2, 4, 6);
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
subplot(2, 4, 7);
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

%saccMore incorr
subplot(2, 4, 8);
hold on;
PlotPSTH(time_axis,TP_ALL_Eff_SaccMore2_Incorrect_Bins_Avg,colorsC(1,:));
PlotPSTH(time_axis,TA_ALL_Eff_SaccMore2_Incorrect_Bins_Avg,colorsC(2,:));
% Labeling
xlabel('Time (s)');
ylabel('Firing Rate (spikes/s)');
title('Efficient search with SaccNum>2 and Incorrect');
legend('TP', 'TA');
grid on;
hold off;


%% Plotting for InEff Correct/Incorrect

time_axis = linspace(-0.6, 1.0, 1600);
colorsC = get_distinguishable_colors(2); % You can also use other colormap functions like parula, jet, etc.
 

figure; % Creates a new figure window
 % Create subplot in a 2x4 grid at the i-th position
%Reject
subplot(2, 4, 1);
hold on;
PlotPSTH(time_axis,TP_ALL_Ineff_Reject_correct_Bins_Avg,colorsC(1,:));
PlotPSTH(time_axis,TA_ALL_Ineff_Reject_correct_Bins_Avg,colorsC(2,:));
% Labeling
xlabel('Time (s)');
ylabel('Firing Rate (spikes/s)');
title('InEfficient search with Reject and Correct');
legend('TP', 'TA');
grid on;
hold off;


%sacc1
subplot(2, 4, 2);
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
subplot(2, 4, 3);
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

%saccMore
subplot(2, 4, 4);
hold on;
PlotPSTH(time_axis,TP_ALL_Ineff_SaccMore2_correct_Bins_Avg,colorsC(1,:));
PlotPSTH(time_axis,TA_ALL_Ineff_SaccMore2_correct_Bins_Avg,colorsC(2,:));
% Labeling
xlabel('Time (s)');
ylabel('Firing Rate (spikes/s)');
title('InEfficient search with SaccNum>2 and Correct');
legend('TP', 'TA');
grid on;
hold off;

%Rej incorr
subplot(2, 4, 5);
hold on;
PlotPSTH(time_axis,TP_ALL_Ineff_Reject_Incorrect_Bins_Avg,colorsC(1,:));
PlotPSTH(time_axis,TA_ALL_Ineff_Reject_Incorrect_Bins_Avg,colorsC(2,:));
% Labeling
xlabel('Time (s)');
ylabel('Firing Rate (spikes/s)');
title('InEfficient search with Reject and Incorrect');
legend('TP', 'TA');
grid on;
hold off;

%sacc1 incorr
subplot(2, 4, 6);
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
subplot(2, 4, 7);
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

%saccMore incorr
subplot(2, 4, 8);
hold on;
PlotPSTH(time_axis,TP_ALL_Ineff_SaccMore2_Incorrect_Bins_Avg,colorsC(1,:));
PlotPSTH(time_axis,TA_ALL_Ineff_SaccMore2_Incorrect_Bins_Avg,colorsC(2,:));
% Labeling
xlabel('Time (s)');
ylabel('Firing Rate (spikes/s)');
title('InEfficient search with SaccNum>2 and Incorrect');
legend('TP', 'TA');
grid on;
hold off;



%% Functions
% Plot PSTH
% Gives mean bins of the data 
function PlotPSTH(time_axis,meanbins,color)

sigma1 = 3; % Adjust sigma value as needed
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





