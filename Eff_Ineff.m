%This code is organized to solve the efficient and inefficient problem

%PSTH plot 
% Get the column names
% Define the column names array
load("SNr_new.mat");



column_names = cell(1, 1600);

for i = 1:1600
    column_names{i} = ['bin', num2str(i)];
end


% Extract data for bins 1 to 1600
bins_1_to_1600_data = table{1, column_names(1:1600)};


% NeuronNumber=20;
% Temp = find(table.iUnit==NeuronNumber);
% Neuron = table(Temp,:);
% 
% Neuronbins = Neuron{:, column_names(1:1600)};
% MeanBins = mean(Neuronbins, 1);


%% TP and TA seperation (2 curves) Mean of neurons


TempTP=find(table.EventValue==4);
TempTA=find(table.EventValue==3);

Eff=find(table.Search_Type==1);
Ineff=find(table.Search_Type==0);

TP_ALL_Eff_temp=intersect(TempTP,Eff);
TP_ALL_Ineff_temp=intersect(TempTP,Ineff);


TA_ALL_Eff_temp=intersect(TempTA,Eff);
TA_ALL_Ineff_temp=intersect(TempTA,Ineff);



TP_ALL_Eff=table(TP_ALL_Eff_temp,:);
TP_ALL_Ineff=table(TP_ALL_Ineff_temp,:);


TA_ALL_Eff=table(TA_ALL_Eff_temp,:);
TA_ALL_Ineff=table(TA_ALL_Ineff_temp,:);





TP_ALL_Eff_Bins=TP_ALL_Eff{:, column_names(1:1600)};
TP_ALL_Ineff_Bins=TP_ALL_Ineff{:, column_names(1:1600)};
TA_ALL_Eff_Bins=TA_ALL_Eff{:, column_names(1:1600)};
TA_ALL_Ineff_Bins=TA_ALL_Ineff{:, column_names(1:1600)};





TA_ALL_Eff_MeanBins=nanmean(TA_ALL_Eff_Bins, 1);
TA_ALL_Ineff_MeanBins=nanmean(TA_ALL_Ineff_Bins, 1);

TP_ALL_Eff_MeanBins=nanmean(TP_ALL_Eff_Bins, 1);
TP_ALL_Ineff_MeanBins=nanmean(TP_ALL_Ineff_Bins, 1);



time_axis = linspace(-0.4, 1.0, 1600);

% Plot firing rates for each condition

figure;
hold on;
colorsC = get_distinguishable_colors(2); % You can also use other colormap functions like parula, jet, etc.
 
%Efficient type of trials

PlotPSTH(time_axis,TP_ALL_Eff_MeanBins,colorsC(1,:));
PlotPSTH(time_axis,TA_ALL_Eff_MeanBins,colorsC(2,:));


% Labeling
xlabel('Time (s)');
ylabel('Firing Rate (spikes/s)');
title('Efficient PSTH Firing Rates for All neurons(n=49)');
legend('TP', 'TA');
grid on;
hold off;



%%
figure;
hold on;

PlotPSTH(time_axis,TP_ALL_Ineff_MeanBins,colorsC(1,:));
PlotPSTH(time_axis,TA_ALL_Ineff_MeanBins,colorsC(2,:));


% Labeling
xlabel('Time (s)');
ylabel('Firing Rate (spikes/s)');
title('Inefficient PSTH Firing Rates for All neurons(n=49)');
legend('TP', 'TA');


grid on;
hold off;




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



