% Checkig the value signals for different conditions SNr

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

%% % TP/TA separation
tableB=table;
TempTP=find(tableB.EventValue==4);
TempTA=find(tableB.EventValue==3);

% Eff/Ineff separation
Eff=find(tableB.Search_Type==1);
Ineff=find(tableB.Search_Type==0);

% Eff=find(tableB.SlopeTP<=20);
% Ineff=find(tableB.SlopeTP>=35);

% Corr/Incorr Separating
Correct=find(tableB.Accuracy==1);
InCorrect=find(tableB.Accuracy==0);

% Sacc Numbers separation
Reject=find(isnan(tableB.SaccNum));
Sacc1=find(tableB.SaccNum==1);
Sacc2=find(tableB.SaccNum==2);
SaccMore2=find(tableB.SaccNum>=2);

% Reject Correct/Incorrect
TP_ALL_Eff_Reject_temp=intersect(intersect(TempTP,Eff),Reject);
TP_ALL_Eff_Reject_temp_correct=intersect(TP_ALL_Eff_Reject_temp,Correct);
TP_ALL_Eff_Reject_temp_Incorrect=intersect(TP_ALL_Eff_Reject_temp,InCorrect);
TP_ALL_Eff_Reject_Incorrect=tableB(TP_ALL_Eff_Reject_temp_Incorrect,:);
TP_ALL_Eff_Reject_correct=tableB(TP_ALL_Eff_Reject_temp_correct,:);
TP_ALL_Eff_Reject_Incorrect_Bins=TP_ALL_Eff_Reject_Incorrect{:, column_names(1:1600)};
TP_ALL_Eff_Reject_correct_Bins=TP_ALL_Eff_Reject_correct{:, column_names(1:1600)};
TP_ALL_Eff_Reject_Incorrect_Bins_Avg=nanmean(TP_ALL_Eff_Reject_Incorrect_Bins,1);%
TP_ALL_Eff_Reject_correct_Bins_Avg=nanmean(TP_ALL_Eff_Reject_correct_Bins,1);%

TP_ALL_Ineff_Reject_temp=intersect(intersect(TempTP,Ineff),Reject);
TP_ALL_Ineff_Reject_temp_correct=intersect(TP_ALL_Ineff_Reject_temp,Correct);
TP_ALL_Ineff_Reject_temp_Incorrect=intersect(TP_ALL_Ineff_Reject_temp,InCorrect);
TP_ALL_Ineff_Reject_Incorrect=tableB(TP_ALL_Ineff_Reject_temp_Incorrect,:);
TP_ALL_Ineff_Reject_correct=tableB(TP_ALL_Ineff_Reject_temp_correct,:);
TP_ALL_Ineff_Reject_Incorrect_Bins=TP_ALL_Ineff_Reject_Incorrect{:, column_names(1:1600)};
TP_ALL_Ineff_Reject_correct_Bins=TP_ALL_Ineff_Reject_correct{:, column_names(1:1600)};
TP_ALL_Ineff_Reject_Incorrect_Bins_Avg=nanmean(TP_ALL_Ineff_Reject_Incorrect_Bins,1);%
TP_ALL_Ineff_Reject_correct_Bins_Avg=nanmean(TP_ALL_Ineff_Reject_correct_Bins,1);%

TA_ALL_Eff_Reject_temp=intersect(intersect(TempTA,Eff),Reject);
TA_ALL_Eff_Reject_temp_correct=intersect(TA_ALL_Eff_Reject_temp,Correct);
TA_ALL_Eff_Reject_temp_Incorrect=intersect(TA_ALL_Eff_Reject_temp,InCorrect);
TA_ALL_Eff_Reject_Incorrect=tableB(TA_ALL_Eff_Reject_temp_Incorrect,:);
TA_ALL_Eff_Reject_correct=tableB(TA_ALL_Eff_Reject_temp_correct,:);
TA_ALL_Eff_Reject_Incorrect_Bins=TA_ALL_Eff_Reject_Incorrect{:, column_names(1:1600)};
TA_ALL_Eff_Reject_correct_Bins=TA_ALL_Eff_Reject_correct{:, column_names(1:1600)};
TA_ALL_Eff_Reject_Incorrect_Bins_Avg=nanmean(TA_ALL_Eff_Reject_Incorrect_Bins,1);%
TA_ALL_Eff_Reject_correct_Bins_Avg=nanmean(TA_ALL_Eff_Reject_correct_Bins,1);%

TA_ALL_Ineff_Reject_temp=intersect(intersect(TempTA,Ineff),Reject);
TA_ALL_Ineff_Reject_temp_correct=intersect(TA_ALL_Ineff_Reject_temp,Correct);
TA_ALL_Ineff_Reject_temp_Incorrect=intersect(TA_ALL_Ineff_Reject_temp,InCorrect);
TA_ALL_Ineff_Reject_Incorrect=tableB(TA_ALL_Ineff_Reject_temp_Incorrect,:);
TA_ALL_Ineff_Reject_correct=tableB(TA_ALL_Ineff_Reject_temp_correct,:);
TA_ALL_Ineff_Reject_Incorrect_Bins=TA_ALL_Ineff_Reject_Incorrect{:, column_names(1:1600)};
TA_ALL_Ineff_Reject_correct_Bins=TA_ALL_Ineff_Reject_correct{:, column_names(1:1600)};
TA_ALL_Ineff_Reject_Incorrect_Bins_Avg=nanmean(TA_ALL_Ineff_Reject_Incorrect_Bins,1);%
TA_ALL_Ineff_Reject_correct_Bins_Avg=nanmean(TA_ALL_Ineff_Reject_correct_Bins,1);%

% 1 Saccade Correct/Incorrect
TP_ALL_Eff_Sacc1_temp=intersect(intersect(TempTP,Eff),Sacc1);
TP_ALL_Eff_Sacc1_temp_correct=intersect(TP_ALL_Eff_Sacc1_temp,Correct);
TP_ALL_Eff_Sacc1_temp_Incorrect=intersect(TP_ALL_Eff_Sacc1_temp,InCorrect);
TP_ALL_Eff_Sacc1_Incorrect=tableB(TP_ALL_Eff_Sacc1_temp_Incorrect,:);
TP_ALL_Eff_Sacc1_correct=tableB(TP_ALL_Eff_Sacc1_temp_correct,:);
TP_ALL_Eff_Sacc1_Incorrect_Bins=TP_ALL_Eff_Sacc1_Incorrect{:, column_names(1:1600)};
TP_ALL_Eff_Sacc1_correct_Bins=TP_ALL_Eff_Sacc1_correct{:, column_names(1:1600)};
TP_ALL_Eff_Sacc1_Incorrect_Bins_Avg=nanmean(TP_ALL_Eff_Sacc1_Incorrect_Bins,1);%
TP_ALL_Eff_Sacc1_correct_Bins_Avg=nanmean(TP_ALL_Eff_Sacc1_correct_Bins,1);%

TP_ALL_Ineff_Sacc1_temp=intersect(intersect(TempTP,Ineff),Sacc1);
TP_ALL_Ineff_Sacc1_temp_correct=intersect(TP_ALL_Ineff_Sacc1_temp,Correct);
TP_ALL_Ineff_Sacc1_temp_Incorrect=intersect(TP_ALL_Ineff_Sacc1_temp,InCorrect);
TP_ALL_Ineff_Sacc1_Incorrect=tableB(TP_ALL_Ineff_Sacc1_temp_Incorrect,:);
TP_ALL_Ineff_Sacc1_correct=tableB(TP_ALL_Ineff_Sacc1_temp_correct,:);
TP_ALL_Ineff_Sacc1_Incorrect_Bins=TP_ALL_Ineff_Sacc1_Incorrect{:, column_names(1:1600)};
TP_ALL_Ineff_Sacc1_correct_Bins=TP_ALL_Ineff_Sacc1_correct{:, column_names(1:1600)};
TP_ALL_Ineff_Sacc1_Incorrect_Bins_Avg=nanmean(TP_ALL_Ineff_Sacc1_Incorrect_Bins,1);%
TP_ALL_Ineff_Sacc1_correct_Bins_Avg=nanmean(TP_ALL_Ineff_Sacc1_correct_Bins,1);%

TA_ALL_Eff_Sacc1_temp=intersect(intersect(TempTA,Eff),Sacc1);
TA_ALL_Eff_Sacc1_temp_correct=intersect(TA_ALL_Eff_Sacc1_temp,Correct);
TA_ALL_Eff_Sacc1_temp_Incorrect=intersect(TA_ALL_Eff_Sacc1_temp,InCorrect);
TA_ALL_Eff_Sacc1_Incorrect=tableB(TA_ALL_Eff_Sacc1_temp_Incorrect,:);
TA_ALL_Eff_Sacc1_correct=tableB(TA_ALL_Eff_Sacc1_temp_correct,:);
TA_ALL_Eff_Sacc1_Incorrect_Bins=TA_ALL_Eff_Sacc1_Incorrect{:, column_names(1:1600)};
TA_ALL_Eff_Sacc1_correct_Bins=TA_ALL_Eff_Sacc1_correct{:, column_names(1:1600)};
TA_ALL_Eff_Sacc1_Incorrect_Bins_Avg=nanmean(TA_ALL_Eff_Sacc1_Incorrect_Bins,1);%
TA_ALL_Eff_Sacc1_correct_Bins_Avg=nanmean(TA_ALL_Eff_Sacc1_correct_Bins,1);%

TA_ALL_Ineff_Sacc1_temp=intersect(intersect(TempTA,Ineff),Sacc1);
TA_ALL_Ineff_Sacc1_temp_correct=intersect(TA_ALL_Ineff_Sacc1_temp,Correct);
TA_ALL_Ineff_Sacc1_temp_Incorrect=intersect(TA_ALL_Ineff_Sacc1_temp,InCorrect);
TA_ALL_Ineff_Sacc1_Incorrect=tableB(TA_ALL_Ineff_Sacc1_temp_Incorrect,:);
TA_ALL_Ineff_Sacc1_correct=tableB(TA_ALL_Ineff_Sacc1_temp_correct,:);
TA_ALL_Ineff_Sacc1_Incorrect_Bins=TA_ALL_Ineff_Sacc1_Incorrect{:, column_names(1:1600)};
TA_ALL_Ineff_Sacc1_correct_Bins=TA_ALL_Ineff_Sacc1_correct{:, column_names(1:1600)};
TA_ALL_Ineff_Sacc1_Incorrect_Bins_Avg=nanmean(TA_ALL_Ineff_Sacc1_Incorrect_Bins,1);%
TA_ALL_Ineff_Sacc1_correct_Bins_Avg=nanmean(TA_ALL_Ineff_Sacc1_correct_Bins,1);%


% 2 Saccade Correct/Incorrect
TP_ALL_Eff_Sacc2_temp=intersect(intersect(TempTP,Eff),Sacc2);
TP_ALL_Eff_Sacc2_temp_correct=intersect(TP_ALL_Eff_Sacc2_temp,Correct);
TP_ALL_Eff_Sacc2_temp_Incorrect=intersect(TP_ALL_Eff_Sacc2_temp,InCorrect);
TP_ALL_Eff_Sacc2_Incorrect=tableB(TP_ALL_Eff_Sacc2_temp_Incorrect,:);
TP_ALL_Eff_Sacc2_correct=tableB(TP_ALL_Eff_Sacc2_temp_correct,:);
TP_ALL_Eff_Sacc2_Incorrect_Bins=TP_ALL_Eff_Sacc2_Incorrect{:, column_names(1:1600)};
TP_ALL_Eff_Sacc2_correct_Bins=TP_ALL_Eff_Sacc2_correct{:, column_names(1:1600)};
TP_ALL_Eff_Sacc2_Incorrect_Bins_Avg=nanmean(TP_ALL_Eff_Sacc2_Incorrect_Bins,1);%
TP_ALL_Eff_Sacc2_correct_Bins_Avg=nanmean(TP_ALL_Eff_Sacc2_correct_Bins,1);%

TP_ALL_Ineff_Sacc2_temp=intersect(intersect(TempTP,Ineff),Sacc2);
TP_ALL_Ineff_Sacc2_temp_correct=intersect(TP_ALL_Ineff_Sacc2_temp,Correct);
TP_ALL_Ineff_Sacc2_temp_Incorrect=intersect(TP_ALL_Ineff_Sacc2_temp,InCorrect);
TP_ALL_Ineff_Sacc2_Incorrect=tableB(TP_ALL_Ineff_Sacc2_temp_Incorrect,:);
TP_ALL_Ineff_Sacc2_correct=tableB(TP_ALL_Ineff_Sacc2_temp_correct,:);
TP_ALL_Ineff_Sacc2_Incorrect_Bins=TP_ALL_Ineff_Sacc2_Incorrect{:, column_names(1:1600)};
TP_ALL_Ineff_Sacc2_correct_Bins=TP_ALL_Ineff_Sacc2_correct{:, column_names(1:1600)};
TP_ALL_Ineff_Sacc2_Incorrect_Bins_Avg=nanmean(TP_ALL_Ineff_Sacc2_Incorrect_Bins,1);%
TP_ALL_Ineff_Sacc2_correct_Bins_Avg=nanmean(TP_ALL_Ineff_Sacc2_correct_Bins,1);%

TA_ALL_Eff_Sacc2_temp=intersect(intersect(TempTA,Eff),Sacc2);
TA_ALL_Eff_Sacc2_temp_correct=intersect(TA_ALL_Eff_Sacc2_temp,Correct);
TA_ALL_Eff_Sacc2_temp_Incorrect=intersect(TA_ALL_Eff_Sacc2_temp,InCorrect);
TA_ALL_Eff_Sacc2_Incorrect=tableB(TA_ALL_Eff_Sacc2_temp_Incorrect,:);
TA_ALL_Eff_Sacc2_correct=tableB(TA_ALL_Eff_Sacc2_temp_correct,:);
TA_ALL_Eff_Sacc2_Incorrect_Bins=TA_ALL_Eff_Sacc2_Incorrect{:, column_names(1:1600)};
TA_ALL_Eff_Sacc2_correct_Bins=TA_ALL_Eff_Sacc2_correct{:, column_names(1:1600)};
TA_ALL_Eff_Sacc2_Incorrect_Bins_Avg=nanmean(TA_ALL_Eff_Sacc2_Incorrect_Bins,1);%
TA_ALL_Eff_Sacc2_correct_Bins_Avg=nanmean(TA_ALL_Eff_Sacc2_correct_Bins,1);%

TA_ALL_Ineff_Sacc2_temp=intersect(intersect(TempTA,Ineff),Sacc2);
TA_ALL_Ineff_Sacc2_temp_correct=intersect(TA_ALL_Ineff_Sacc2_temp,Correct);
TA_ALL_Ineff_Sacc2_temp_Incorrect=intersect(TA_ALL_Ineff_Sacc2_temp,InCorrect);
TA_ALL_Ineff_Sacc2_Incorrect=tableB(TA_ALL_Ineff_Sacc2_temp_Incorrect,:);
TA_ALL_Ineff_Sacc2_correct=tableB(TA_ALL_Ineff_Sacc2_temp_correct,:);
TA_ALL_Ineff_Sacc2_Incorrect_Bins=TA_ALL_Ineff_Sacc2_Incorrect{:, column_names(1:1600)};
TA_ALL_Ineff_Sacc2_correct_Bins=TA_ALL_Ineff_Sacc2_correct{:, column_names(1:1600)};
TA_ALL_Ineff_Sacc2_Incorrect_Bins_Avg=nanmean(TA_ALL_Ineff_Sacc2_Incorrect_Bins,1);%
TA_ALL_Ineff_Sacc2_correct_Bins_Avg=nanmean(TA_ALL_Ineff_Sacc2_correct_Bins,1);%



% More than 2 Saccade Correct/Incorrect
TP_ALL_Eff_SaccMore2_temp=intersect(intersect(TempTP,Eff),SaccMore2);
TP_ALL_Eff_SaccMore2_temp_correct=intersect(TP_ALL_Eff_SaccMore2_temp,Correct);
TP_ALL_Eff_SaccMore2_temp_Incorrect=intersect(TP_ALL_Eff_SaccMore2_temp,InCorrect);
TP_ALL_Eff_SaccMore2_Incorrect=tableB(TP_ALL_Eff_SaccMore2_temp_Incorrect,:);
TP_ALL_Eff_SaccMore2_correct=tableB(TP_ALL_Eff_SaccMore2_temp_correct,:);
TP_ALL_Eff_SaccMore2_Incorrect_Bins=TP_ALL_Eff_SaccMore2_Incorrect{:, column_names(1:1600)};
TP_ALL_Eff_SaccMore2_correct_Bins=TP_ALL_Eff_SaccMore2_correct{:, column_names(1:1600)};
TP_ALL_Eff_SaccMore2_Incorrect_Bins_Avg=nanmean(TP_ALL_Eff_SaccMore2_Incorrect_Bins,1);%
TP_ALL_Eff_SaccMore2_correct_Bins_Avg=nanmean(TP_ALL_Eff_SaccMore2_correct_Bins,1);%

TP_ALL_Ineff_SaccMore2_temp=intersect(intersect(TempTP,Ineff),SaccMore2);
TP_ALL_Ineff_SaccMore2_temp_correct=intersect(TP_ALL_Ineff_SaccMore2_temp,Correct);
TP_ALL_Ineff_SaccMore2_temp_Incorrect=intersect(TP_ALL_Ineff_SaccMore2_temp,InCorrect);
TP_ALL_Ineff_SaccMore2_Incorrect=tableB(TP_ALL_Ineff_SaccMore2_temp_Incorrect,:);
TP_ALL_Ineff_SaccMore2_correct=tableB(TP_ALL_Ineff_SaccMore2_temp_correct,:);
TP_ALL_Ineff_SaccMore2_Incorrect_Bins=TP_ALL_Ineff_SaccMore2_Incorrect{:, column_names(1:1600)};
TP_ALL_Ineff_SaccMore2_correct_Bins=TP_ALL_Ineff_SaccMore2_correct{:, column_names(1:1600)};
TP_ALL_Ineff_SaccMore2_Incorrect_Bins_Avg=nanmean(TP_ALL_Ineff_SaccMore2_Incorrect_Bins,1);%
TP_ALL_Ineff_SaccMore2_correct_Bins_Avg=nanmean(TP_ALL_Ineff_SaccMore2_correct_Bins,1);%

TA_ALL_Eff_SaccMore2_temp=intersect(intersect(TempTA,Eff),SaccMore2);
TA_ALL_Eff_SaccMore2_temp_correct=intersect(TA_ALL_Eff_SaccMore2_temp,Correct);
TA_ALL_Eff_SaccMore2_temp_Incorrect=intersect(TA_ALL_Eff_SaccMore2_temp,InCorrect);
TA_ALL_Eff_SaccMore2_Incorrect=tableB(TA_ALL_Eff_SaccMore2_temp_Incorrect,:);
TA_ALL_Eff_SaccMore2_correct=tableB(TA_ALL_Eff_SaccMore2_temp_correct,:);
TA_ALL_Eff_SaccMore2_Incorrect_Bins=TA_ALL_Eff_SaccMore2_Incorrect{:, column_names(1:1600)};
TA_ALL_Eff_SaccMore2_correct_Bins=TA_ALL_Eff_SaccMore2_correct{:, column_names(1:1600)};
TA_ALL_Eff_SaccMore2_Incorrect_Bins_Avg=nanmean(TA_ALL_Eff_SaccMore2_Incorrect_Bins,1);%
TA_ALL_Eff_SaccMore2_correct_Bins_Avg=nanmean(TA_ALL_Eff_SaccMore2_correct_Bins,1);%

TA_ALL_Ineff_SaccMore2_temp=intersect(intersect(TempTA,Ineff),SaccMore2);
TA_ALL_Ineff_SaccMore2_temp_correct=intersect(TA_ALL_Ineff_SaccMore2_temp,Correct);
TA_ALL_Ineff_SaccMore2_temp_Incorrect=intersect(TA_ALL_Ineff_SaccMore2_temp,InCorrect);
TA_ALL_Ineff_SaccMore2_Incorrect=tableB(TA_ALL_Ineff_SaccMore2_temp_Incorrect,:);
TA_ALL_Ineff_SaccMore2_correct=tableB(TA_ALL_Ineff_SaccMore2_temp_correct,:);
TA_ALL_Ineff_SaccMore2_Incorrect_Bins=TA_ALL_Ineff_SaccMore2_Incorrect{:, column_names(1:1600)};
TA_ALL_Ineff_SaccMore2_correct_Bins=TA_ALL_Ineff_SaccMore2_correct{:, column_names(1:1600)};
TA_ALL_Ineff_SaccMore2_Incorrect_Bins_Avg=nanmean(TA_ALL_Ineff_SaccMore2_Incorrect_Bins,1);%
TA_ALL_Ineff_SaccMore2_correct_Bins_Avg=nanmean(TA_ALL_Ineff_SaccMore2_correct_Bins,1);%

%% 



time_axis = linspace(-0.6, 1.0, 1600);
colorsC = get_distinguishable_colors(2); % You can also use other colormap functions like parula, jet, etc.
 

figure; % Creates a new figure window
 % Create subplot in a 2x4 grid at the i-th position


%sacc2
subplot(2, 2, 1);
hold on;
ValueSignal_Eff_Sacc2_correct=-(TA_ALL_Eff_Sacc2_correct_Bins_Avg-TP_ALL_Eff_Sacc2_correct_Bins_Avg);
ValueSignal_Eff_Sacc2_Incorrect=-(TA_ALL_Eff_Sacc2_Incorrect_Bins_Avg-TP_ALL_Eff_Sacc2_Incorrect_Bins_Avg);
PlotPSTH(time_axis,ValueSignal_Eff_Sacc2_correct,colorsC(1,:));
PlotPSTH(time_axis,ValueSignal_Eff_Sacc2_Incorrect,colorsC(2,:));
% Labeling
xlabel('Time (s)');
ylabel('Firing Rate (spikes/s)');
title('Value signal for Incorrect/correct Efficient Sacc2');
legend('Correct', 'Incorrect');
ax = gca;
ax.FontSize = 19;
grid on;
hold off;



%saccMore
subplot(2, 2, 2);
hold on;
ValueSignal_Eff_SaccMore2_correct=-(TA_ALL_Eff_SaccMore2_correct_Bins_Avg-TP_ALL_Eff_SaccMore2_correct_Bins_Avg);
ValueSignal_Eff_SaccMore2_Incorrect=-(TA_ALL_Eff_SaccMore2_Incorrect_Bins_Avg-TP_ALL_Eff_SaccMore2_Incorrect_Bins_Avg);
PlotPSTH(time_axis,ValueSignal_Eff_SaccMore2_correct,colorsC(1,:));
PlotPSTH(time_axis,ValueSignal_Eff_SaccMore2_Incorrect,colorsC(2,:));
% Labeling
xlabel('Time (s)');
ylabel('Firing Rate (spikes/s)');
title('Value signal for Incorrect/correct Efficient Sacc>=2');
legend('Correct', 'Incorrect');
ax = gca;
ax.FontSize = 19;
grid on;
hold off;



%sacc2
subplot(2, 2, 3);
hold on;
ValueSignal_Ineeficient_Sacc2_correct=-(TA_ALL_Ineff_Sacc2_correct_Bins_Avg-TP_ALL_Ineff_Sacc2_correct_Bins_Avg);
ValueSignal_Ineeficient_Sacc2_Incorrect=-(TA_ALL_Ineff_Sacc2_Incorrect_Bins_Avg-TP_ALL_Ineff_Sacc2_Incorrect_Bins_Avg);
PlotPSTH(time_axis,ValueSignal_Ineeficient_Sacc2_correct,colorsC(1,:));
PlotPSTH(time_axis,ValueSignal_Ineeficient_Sacc2_Incorrect,colorsC(2,:));
% Labeling
xlabel('Time (s)');
ylabel('Firing Rate (spikes/s)');
title('Value signal for Incorrect/correct InEfficient Sacc2');
legend('Correct', 'Incorrect');
ax = gca;
ax.FontSize = 19;
grid on;
hold off;




%saccMore
subplot(2, 2, 4);
hold on;
ValueSignal_Ineff_SaccMore2_Correct=-(TA_ALL_Ineff_SaccMore2_correct_Bins_Avg-TP_ALL_Ineff_SaccMore2_correct_Bins_Avg);
ValueSignal_Ineff_SaccMore2_InCorrect=-(TA_ALL_Ineff_SaccMore2_Incorrect_Bins_Avg-TP_ALL_Ineff_SaccMore2_Incorrect_Bins_Avg);
PlotPSTH(time_axis,ValueSignal_Ineff_SaccMore2_Correct,colorsC(1,:));
PlotPSTH(time_axis,ValueSignal_Ineff_SaccMore2_InCorrect,colorsC(2,:));
% Labeling
xlabel('Time (s)');
ylabel('Firing Rate (spikes/s)');
title('Value signal for Incorrect/correct InEfficient Sacc>=2');
legend('Correct', 'Incorrect');
ax = gca;
ax.FontSize = 19;
grid on;
hold off;


%% All of Correcct/Incorrect

TempTP=find(tableB.EventValue==4);
TempTA=find(tableB.EventValue==3);

% Eff/Ineff separation
%Eff=find(tableB.Search_Type==1);
%Ineff=find(tableB.Search_Type==0);

 Eff=find(tableB.SlopeTP<=20);
 Ineff=find(tableB.SlopeTP>=35);

% Corr/Incorr Separating
Correct=find(tableB.Accuracy==1);
InCorrect=find(tableB.Accuracy==0);

%Eff
%Corr

TP_Eff_Correct_temp=intersect(intersect(TempTP,Eff),Correct);
TP_Eff_Correct=tableB(TP_Eff_Correct_temp,:);
TP_Eff_Correct_Bins=TP_Eff_Correct{:, column_names(1:1600)};
TP_Eff_Correct_Bins_avrg=nanmean(TP_Eff_Correct_Bins,1);

TA_Eff_only=intersect(TempTA,Eff);
TA_Eff_Correct_temp=intersect(intersect(TempTA,Eff),Correct);
TA_Eff_Correct=tableB(TA_Eff_Correct_temp,:);
TA_Eff_Correct_Bins=TA_Eff_Correct{:, column_names(1:1600)};
TA_Eff_Correct_Bins_avrg=nanmean(TA_Eff_Correct_Bins,1);

%Value Signal
ValueSIgnal_Eff_correct=-(TA_Eff_Correct_Bins_avrg-TP_Eff_Correct_Bins_avrg);

%Incorr
TP_Eff_InCorrect_temp=intersect(intersect(TempTP,Eff),InCorrect);
TP_Eff_InCorrect=tableB(TP_Eff_InCorrect_temp,:);
TP_Eff_InCorrect_Bins=TP_Eff_InCorrect{:, column_names(1:1600)};
TP_Eff_InCorrect_Bins_avrg=nanmean(TP_Eff_InCorrect_Bins,1);

TA_Eff_InCorrect_temp=intersect(intersect(TempTA,Eff),InCorrect);
TA_Eff_InCorrect=tableB(TA_Eff_InCorrect_temp,:);
TA_Eff_InCorrect_Bins=TA_Eff_InCorrect{:, column_names(1:1600)};
TA_Eff_InCorrect_Bins_avrg=nanmean(TA_Eff_InCorrect_Bins,1);

%Value Signal
ValueSIgnal_Eff_Incorrect=-(TA_Eff_InCorrect_Bins_avrg-TP_Eff_InCorrect_Bins_avrg);

%InEFF
%Corr
TP_Ineff=intersect(TempTP,Ineff);

TP_InEff_Correct_temp=intersect(intersect(TempTP,Ineff),Correct);
TP_InEff_Correct=tableB(TP_InEff_Correct_temp,:);
TP_InEff_Correct_Bins=TP_InEff_Correct{:, column_names(1:1600)};
TP_InEff_Correct_Bins_avrg=nanmean(TP_InEff_Correct_Bins,1);

TA_InEff=intersect(TempTA,Ineff);
TA_InEff_Correct_temp=intersect(intersect(TempTA,Ineff),Correct);
TA_InEff_Correct=tableB(TA_InEff_Correct_temp,:);
TA_InEff_Correct_Bins=TA_InEff_Correct{:, column_names(1:1600)};
TA_InEff_Correct_Bins_avrg=nanmean(TA_InEff_Correct_Bins,1);

%Value Signal
ValueSignal_InEff_Correct=-(TA_InEff_Correct_Bins_avrg-TP_InEff_Correct_Bins_avrg);




%Incorr
TP_InEff_InCorrect_temp=intersect(intersect(TempTP,Ineff),InCorrect);
TP_InEff_InCorrect=tableB(TP_InEff_InCorrect_temp,:);
TP_InEff_InCorrect_Bins=TP_InEff_InCorrect{:, column_names(1:1600)};
TP_InEff_InCorrect_Bins_avrg=nanmean(TP_InEff_InCorrect_Bins,1);

TA_InEff_InCorrect_temp=intersect(intersect(TempTA,Ineff),InCorrect);
TA_InEff_InCorrect=tableB(TA_InEff_InCorrect_temp,:);
TA_InEff_InCorrect_Bins=TA_InEff_InCorrect{:, column_names(1:1600)};
TA_InEff_InCorrect_Bins_avrg=nanmean(TA_InEff_InCorrect_Bins,1);

%Value Signal
ValueSignal_InEff_InCorrect=-(TA_InEff_InCorrect_Bins_avrg-TP_InEff_InCorrect_Bins_avrg);



time_axis = linspace(-0.6, 1.0, 1600);
colorsC = get_distinguishable_colors(2); % You can also use other colormap functions like parula, jet, etc.
 

figure; % Creates a new figure window
 % Create subplot in a 2x4 grid at the i-th position


%sacc2
subplot(2, 1, 1);
hold on;

PlotPSTH(time_axis,ValueSIgnal_Eff_correct,colorsC(1,:));
PlotPSTH(time_axis,ValueSIgnal_Eff_Incorrect,colorsC(2,:));
% Labeling
xlabel('Time (s)');
ylabel('Firing Rate (spikes/s)');
title('Value signal for Incorrect/correct Efficient');
legend('Correct', 'Incorrect');
ax = gca;
ax.FontSize = 19;
grid on;
hold off;



%saccMore
subplot(2, 1, 2);
hold on;
PlotPSTH(time_axis,ValueSignal_InEff_Correct,colorsC(1,:));
PlotPSTH(time_axis,ValueSignal_InEff_InCorrect,colorsC(2,:));
% Labeling
xlabel('Time (s)');
ylabel('Firing Rate (spikes/s)');
title('Value signal for Incorrect/correct InEfficient ');
legend('Correct', 'Incorrect');
ax = gca;
ax.FontSize = 19;
grid on;
hold off;


%% value sig eff ineff



TP_Eff_Only_temp=intersect(TempTP,Eff);
TP_Eff_Only=tableB(TP_Eff_Only_temp,:);
TP_Eff_Only_Bins=TP_Eff_Only{:, column_names(1:1600)};
TP_Eff_Only_Bins_avrg=nanmean(TP_Eff_Only_Bins,1);

TA_Eff_Only_temp=intersect(TempTA,Eff);
TA_Eff_Only=tableB(TA_Eff_Only_temp,:);
TA_Eff_Only_Bins=TA_Eff_Only{:, column_names(1:1600)};
TA_Eff_Only_Bins_avrg=nanmean(TA_Eff_Only_Bins,1);

valueSignal_Eff=TP_Eff_Only_Bins_avrg-TA_Eff_Only_Bins_avrg;



TP_Ineff_Only_temp=intersect(TempTP,Ineff);
TP_Ineff_Only=tableB(TP_Ineff_Only_temp,:);
TP_Ineff_Only_Bins=TP_Ineff_Only{:, column_names(1:1600)};
TP_Ineff_Only_Bins_avrg=nanmean(TP_Ineff_Only_Bins,1);

TA_Ineff_Only_temp=intersect(TempTA,Ineff);
TA_Ineff_Only=tableB(TA_Ineff_Only_temp,:);
TA_Ineff_Only_Bins=TA_Ineff_Only{:, column_names(1:1600)};
TA_Ineff_Only_Bins_avrg=nanmean(TA_Ineff_Only_Bins,1);

ValueSignal_Ineff=TP_Ineff_Only_Bins_avrg-TA_Ineff_Only_Bins_avrg;

time_axis = linspace(-0.6, 1.0, 1600);
colorsC = get_distinguishable_colors(2); % You can also use other colormap functions like parula, jet, etc.
 

figure; % Creates a new figure window
 % Create subplot in a 2x4 grid at the i-th position


%sacc2
subplot(2, 1, 1);
hold on;

PlotPSTH(time_axis,valueSignal_Eff,colorsC(1,:));
% Labeling
xlabel('Time (s)');
ylabel('Firing Rate (spikes/s)');
title('Value signal for  Efficient');
ax = gca;
ax.FontSize = 19;
grid on;
hold off;



%saccMore
subplot(2, 1, 2);
hold on;
PlotPSTH(time_axis,ValueSignal_Ineff,colorsC(1,:));
% Labeling
xlabel('Time (s)');
ylabel('Firing Rate (spikes/s)');
title('Value signal for  InEfficient ');
legend('Correct', 'Incorrect');
ax = gca;
ax.FontSize = 19;
grid on;
hold off;



%% Functions
% Plot PSTH
% Gives mean bins of the data 
function PlotPSTH(time_axis,meanbins,color)

sigma1 = 15; % Adjust sigma value as needed
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





