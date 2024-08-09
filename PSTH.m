%PSTH plot 
% Get the column names
% Define the column names array
load("SNr.mat")
load('LFP_Spike_SNr.mat');

load("vlPFC.mat")
Temp2 = find(T.iUnit==22);
Neuron2 = T(Temp2,:);



column_names = cell(1, 1600);

for i = 1:1600
    column_names{i} = ['bin', num2str(i)];
end


% Extract data for bins 1 to 1600
bins_1_to_1600_data = table{1, column_names(1:1600)};


NeuronNumber=20;
Temp = find(table.iUnit==NeuronNumber);
Neuron = table(Temp,:);

Neuronbins = Neuron{:, column_names(1:1600)};
MeanBins = mean(Neuronbins, 1);





%% Seperation for Neuron #N=20
%Neuron=table;
Neuron=tableB;

TP = find(Neuron.EventValue==4);
TA = find(Neuron.EventValue==3);

DspSize3 = find(Neuron.DispSize==3);
DspSize5 = find(Neuron.DispSize==5);
DspSize7 = find(Neuron.DispSize==7);
DspSize9 = find(Neuron.DispSize==9);

%efficient=find(Neuron.Efficient=='1');
%inefficient=find(Neuron2.Efficient=='0');


TP_Dspsz3=intersect(TP,DspSize3);
Neuron_TP_Dspsz3=Neuron(TP_Dspsz3,:);
Neuron_TP_Dspsz3_bins = Neuron_TP_Dspsz3{:, column_names(1:1600)};
Neuron_TP_Dspsz3_bins_mean=nanmean(Neuron_TP_Dspsz3_bins,1);

TP_Dspsz5=intersect(TP,DspSize5);
Neuron_TP_Dspsz5=Neuron(TP_Dspsz5,:);
Neuron_TP_Dspsz5_bins = Neuron_TP_Dspsz5{:, column_names(1:1600)};
Neuron_TP_Dspsz5_bins_mean=nanmean(Neuron_TP_Dspsz5_bins,1);


TP_Dspsz7=intersect(TP,DspSize7);
Neuron_TP_Dspsz7=Neuron(TP_Dspsz7,:);
Neuron_TP_Dspsz7_bins = Neuron_TP_Dspsz7{:, column_names(1:1600)};
Neuron_TP_Dspsz7_bins_mean=nanmean(Neuron_TP_Dspsz7_bins,1);


TP_Dspsz9=intersect(TP,DspSize9);
Neuron_TP_Dspsz9=Neuron(TP_Dspsz9,:);
Neuron_TP_Dspsz9_bins = Neuron_TP_Dspsz9{:, column_names(1:1600)};
Neuron_TP_Dspsz9_bins_mean=nanmean(Neuron_TP_Dspsz9_bins,1);


TA_Dspsz3=intersect(TA,DspSize3);
Neuron_TA_Dspsz3=Neuron(TA_Dspsz3,:);
Neuron_TA_Dspsz3_bins = Neuron_TA_Dspsz3{:, column_names(1:1600)};
Neuron_TA_Dspsz3_bins_mean=nanmean(Neuron_TA_Dspsz3_bins,1);


TA_Dspsz5=intersect(TA,DspSize5);
Neuron_TA_Dspsz5=Neuron(TA_Dspsz5,:);
Neuron_TA_Dspsz5_bins = Neuron_TA_Dspsz5{:, column_names(1:1600)};
Neuron_TA_Dspsz5_bins_mean=nanmean(Neuron_TA_Dspsz5_bins,1);


TA_Dspsz7=intersect(TA,DspSize7);
Neuron_TA_Dspsz7=Neuron(TA_Dspsz7,:);
Neuron_TA_Dspsz7_bins = Neuron_TA_Dspsz7{:, column_names(1:1600)};
Neuron_TA_Dspsz7_bins_mean=nanmean(Neuron_TA_Dspsz7_bins,1);


TA_Dspsz9=intersect(TA,DspSize9);
Neuron_TA_Dspsz9=Neuron(TA_Dspsz9,:);
Neuron_TA_Dspsz9_bins = Neuron_TA_Dspsz9{:, column_names(1:1600)};
Neuron_TA_Dspsz9_bins_mean=nanmean(Neuron_TA_Dspsz9_bins,1);




%% Plot the seperated curves of previous section
% Assuming firing_rates is a matrix of size 8x1600 containing firing rates for 8 conditions

% Define time axis (assuming your data spans 1 second)
time_axis = linspace(0, 1.6, 1600);

% Plot firing rates for each condition
figure;
hold on;
colors = get_distinguishable_colors(8); % You can also use other colormap functions like parula, jet, etc.
 

%TODO
PlotPSTH(time_axis,Neuron_TP_Dspsz3_bins_mean,colors(1,:));
PlotPSTH(time_axis,Neuron_TP_Dspsz5_bins_mean,colors(2,:));
PlotPSTH(time_axis,Neuron_TP_Dspsz7_bins_mean,colors(3,:));
PlotPSTH(time_axis,Neuron_TP_Dspsz9_bins_mean,colors(4,:));
PlotPSTH(time_axis,Neuron_TA_Dspsz3_bins_mean,colors(5,:));
PlotPSTH(time_axis,Neuron_TA_Dspsz5_bins_mean,colors(6,:));
PlotPSTH(time_axis,Neuron_TA_Dspsz7_bins_mean,colors(7,:));
PlotPSTH(time_axis,Neuron_TA_Dspsz9_bins_mean,colors(8,:));


% Labeling
xlabel('Time (s)');
ylabel('Firing Rate (spikes/s)');
title('PSTH Firing Rates Across Conditions');
legend('TP DSP3', 'TP DSP5', 'TP DSP7', 'TP DSP9', ...
       'TA DSP3', 'TA DSP5', 'TA DSP7', 'TA DSP9');
grid on;
hold off;



%% TP and TA plot for one neuron #N=20

TPOnly = find(Neuron.EventValue==4);
TAOnly = find(Neuron.EventValue==3);



Neuron_TPOnly=Neuron(TPOnly,:);
Neuron_TAOnly=Neuron(TAOnly,:);

Neuron_TPOnly_bins = Neuron_TPOnly{:, column_names(1:1600)};
Neuron_TAOnly_bins = Neuron_TAOnly{:, column_names(1:1600)};

Neuron_TPOnly_bins_mean=nanmean(Neuron_TPOnly_bins,1);
Neuron_TAOnly_bins_mean=nanmean(Neuron_TAOnly_bins,1);



% Plot firing rates for each condition
figure;
hold on;
colorsB = get_distinguishable_colors(2); % You can also use other colormap functions like parula, jet, etc.
 


PlotPSTH(time_axis,Neuron_TPOnly_bins_mean,colorsB(1,:));
PlotPSTH(time_axis,Neuron_TAOnly_bins_mean,colorsB(2,:));

% Labeling
xlabel('Time (s)');
ylabel('Firing Rate (spikes/s)');
title('PSTH Firing Rates Across Conditions');
legend('TP', 'TA');
grid on;
hold off;

%% finding the number of ALL neurons that have PSTH

Boolean=true;
Count=0;
NeuronsNumber=0;
while(Boolean)
    Count=Count+1;
    if ismember(Count, table.iUnit(:,:))
        NeuronsNumber=NeuronsNumber+1;
        Temporary=find(table.iUnit==Count);
        CurrentNeuron=table(Temporary,:);
        if isnan(CurrentNeuron{1,'bin1'})
            disp(NeuronsNumber);
            Boolean=false;
        else
            Boolean=true;
        end
    end
end

%N neurons are 38

%% TP and TA seperation (2 curves) Mean of neurons


TempTP=find(table.EventValue==4);
TempTA=find(table.EventValue==3);



TP_ALL=table(TempTP,:);
TA_ALL=table(TempTA,:);


TP_ALL_Bins=TP_ALL{:, column_names(1:1600)};
TA_ALL_Bins=TA_ALL{:, column_names(1:1600)};

TP_ALL_MeanBins=nanmean(TP_ALL_Bins, 1);
TA_ALL_MeanBins=nanmean(TA_ALL_Bins, 1);

time_axis = linspace(0, 1.6, 1600);

% Plot firing rates for each condition
figure;
hold on;
colorsC = get_distinguishable_colors(2); % You can also use other colormap functions like parula, jet, etc.
 


PlotPSTH(time_axis,TP_ALL_MeanBins,colorsC(1,:));
PlotPSTH(time_axis,TA_ALL_MeanBins,colorsC(2,:));

% Labeling
xlabel('Time (s)');
ylabel('Firing Rate (spikes/s)');
title('PSTH Firing Rates for All neurons(n=38)');
legend('TP', 'TA');
grid on;
hold off;

%% TP/TA and also DSP size seperation(8 curves) for mean of neurons

TempDspSize3_ALL = find(table.DispSize==3);
TempDspSize5_ALL = find(table.DispSize==5);
TempDspSize7_ALL = find(table.DispSize==7);
TempDspSize9_ALL = find(table.DispSize==9);


%DSP size 3
AllTP_Dspsz3=intersect(TempTP,TempDspSize3_ALL);
AllNeuron_TP_Dspsz3=table(AllTP_Dspsz3,:);
AllNeuron_TP_Dspsz3_bins = AllNeuron_TP_Dspsz3{:, column_names(1:1600)};
AllNeuron_TP_Dspsz3_bins_mean=nanmean(AllNeuron_TP_Dspsz3_bins,1);

%DSP size 5
AllTP_Dspsz5=intersect(TempTP,TempDspSize5_ALL);
AllNeuron_TP_Dspsz5=table(AllTP_Dspsz5,:);
AllNeuron_TP_Dspsz5_bins = AllNeuron_TP_Dspsz5{:, column_names(1:1600)};
AllNeuron_TP_Dspsz5_bins_mean=nanmean(AllNeuron_TP_Dspsz5_bins,1);


%DSP size 7
AllTP_Dspsz7=intersect(TempTP,TempDspSize7_ALL);
AllNeuron_TP_Dspsz7=table(AllTP_Dspsz7,:);
AllNeuron_TP_Dspsz7_bins = AllNeuron_TP_Dspsz7{:, column_names(1:1600)};
AllNeuron_TP_Dspsz7_bins_mean=nanmean(AllNeuron_TP_Dspsz7_bins,1);

%DSP size 9
AllTP_Dspsz9=intersect(TempTP,TempDspSize9_ALL);
AllNeuron_TP_Dspsz9=table(AllTP_Dspsz9,:);
AllNeuron_TP_Dspsz9_bins = AllNeuron_TP_Dspsz9{:, column_names(1:1600)};
AllNeuron_TP_Dspsz9_bins_mean=nanmean(AllNeuron_TP_Dspsz9_bins,1);



%DSP size 3
AllTA_Dspsz3=intersect(TempTA,TempDspSize3_ALL);
AllNeuron_TA_Dspsz3=table(AllTA_Dspsz3,:);
AllNeuron_TA_Dspsz3_bins = AllNeuron_TA_Dspsz3{:, column_names(1:1600)};
AllNeuron_TA_Dspsz3_bins_mean=nanmean(AllNeuron_TA_Dspsz3_bins,1);

%DSP size 5
AllTA_Dspsz5=intersect(TempTA,TempDspSize5_ALL);
AllNeuron_TA_Dspsz5=table(AllTA_Dspsz5,:);
AllNeuron_TA_Dspsz5_bins = AllNeuron_TA_Dspsz5{:, column_names(1:1600)};
AllNeuron_TA_Dspsz5_bins_mean=nanmean(AllNeuron_TA_Dspsz5_bins,1);


%DSP size 7
AllTA_Dspsz7=intersect(TempTA,TempDspSize5_ALL);
AllNeuron_TA_Dspsz7=table(AllTA_Dspsz7,:);
AllNeuron_TA_Dspsz7_bins = AllNeuron_TA_Dspsz7{:, column_names(1:1600)};
AllNeuron_TA_Dspsz7_bins_mean=nanmean(AllNeuron_TA_Dspsz7_bins,1);


%DSP size 9
AllTA_Dspsz9=intersect(TempTA,TempDspSize9_ALL);
AllNeuron_TA_Dspsz9=table(AllTA_Dspsz9,:);
AllNeuron_TA_Dspsz9_bins = AllNeuron_TA_Dspsz9{:, column_names(1:1600)};
AllNeuron_TA_Dspsz9_bins_mean=nanmean(AllNeuron_TA_Dspsz9_bins,1);



% Plot firing rates for each condition
figure;
hold on;
colors = get_distinguishable_colors(8); % You can also use other colormap functions like parula, jet, etc.
 


%TODO
PlotPSTH(time_axis,AllNeuron_TP_Dspsz3_bins_mean,colors(1,:));
PlotPSTH(time_axis,AllNeuron_TP_Dspsz5_bins_mean,colors(2,:));
PlotPSTH(time_axis,AllNeuron_TP_Dspsz7_bins_mean,colors(3,:));
PlotPSTH(time_axis,AllNeuron_TP_Dspsz9_bins_mean,colors(4,:));
PlotPSTH(time_axis,AllNeuron_TA_Dspsz3_bins_mean,colors(5,:));
PlotPSTH(time_axis,AllNeuron_TA_Dspsz5_bins_mean,colors(6,:));
PlotPSTH(time_axis,AllNeuron_TA_Dspsz7_bins_mean,colors(7,:));
PlotPSTH(time_axis,AllNeuron_TA_Dspsz9_bins_mean,colors(8,:));


% Labeling
xlabel('Time (s)');
ylabel('Firing Rate (spikes/s)');
title('PSTH Firing Rates Across Conditions For average of all neurons');
legend('TP DSP3', 'TP DSP5', 'TP DSP7', 'TP DSP9', ...
       'TA DSP3', 'TA DSP5', 'TA DSP7', 'TA DSP9');
grid on;
hold off;




%% ---------------------------------------------------------------------------------------------------------
% PLot


% Define time axis (assuming your data spans 1 second)
time_axis = linspace(0, 1, 1600);
psth_data=MeanBins;
% Plot PSTH as a continuous line
plot(time_axis, psth_data, 'LineWidth', 2);

% Labeling
xlabel('Time (s)');
ylabel('Spike Count');
title('Peri-Stimulus Time Histogram (PSTH)');


%% Smoothed PSTH


%Moving Average Filter
window_size = 5; % Adjust window size as needed
smoothed_data1 = movmean(MeanBins, window_size);

%Gaussian Filter
sigma = 2; % Adjust sigma value as needed
smoothed_data2 = imgaussfilt(MeanBins, sigma);



% Define subplot layout: 1 row, 2 columns
subplot(1, 2, 1); % First subplot

% Plot PSTH as a continuous line
plot(time_axis, smoothed_data1, 'LineWidth', 2);

% Labeling
xlabel('Time (s)');
ylabel('Spike Count');
title('PSTH 1');

subplot(1, 2, 2); % Second subplot

% Plot PSTH as a continuous line
plot(time_axis, smoothed_data2, 'LineWidth', 2);

% Labeling
xlabel('Time (s)');
ylabel('Spike Count');
title('PSTH 2');


%% Neuron 22

NeuronNumber2=22;
Temp2 = find(table.iUnit==NeuronNumber);
Neuron2 = table(Temp2,:);

Neuronbins2 = Neuron2{:, column_names(1:1600)};
MeanBins2 = mean(Neuronbins2, 1);


%Moving Average Filter
window_size = 5; % Adjust window size as needed
smoothed_data1p = movmean(MeanBins2, window_size);

%Gaussian Filter
sigma = 2; % Adjust sigma value as needed
smoothed_data2p = imgaussfilt(MeanBins2, sigma);



% Define subplot layout: 1 row, 2 columns
subplot(1, 2, 1); % First subplot

% Plot PSTH as a continuous line
plot(time_axis, smoothed_data1p, 'LineWidth', 2);

% Labeling
xlabel('Time (s)');
ylabel('Spike Count');
title('PSTH 1');

subplot(1, 2, 2); % Second subplot

% Plot PSTH as a continuous line
plot(time_axis, smoothed_data2p, 'LineWidth', 2);

% Labeling
xlabel('Time (s)');
ylabel('Spike Count');
title('PSTH 2');

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

