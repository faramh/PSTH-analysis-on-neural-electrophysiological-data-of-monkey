clc
clear
close all

load('LFP_Spike_SNr.mat');

Temp = find(T.iUnit==2);
Neuron = T(Temp,:);
clear T
%%
DS3 = find(Neuron.DispSize==7);
Eff = find(Neuron.Efficient=='1');
TP = find(Neuron.EventValue==3);
TA = find(Neuron.EventValue==4);

DS3EffTP = intersect(intersect(DS3,Eff),TP);
DS3EffTA = intersect(intersect(DS3,Eff),TA);

for i=1:length(DS3EffTA)
    str.trial(i).times =  Neuron.SpikeTimes{DS3EffTA(i),1};
end
%%
hold all;
for trialCount = 1:size(str.trial,2)
    spikePos = str.trial(trialCount).times  
    for spikeCount = 1:length(spikePos)
        plot([spikePos(spikeCount) spikePos(spikeCount)], ...
            [trialCount-0.4 trialCount+0.4], 'k');
    end
end
xlim([-200,1400])
ylim([0,size(str.trial,2)])