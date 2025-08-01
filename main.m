% This script demonstrates how the model learns optimized interaction weights 
% across a given set of acoustic features (envelope, pitch, low-frequency 
% spectrogram, high-frequency spectrogram, scale, and rate).
%
% The demonstration uses half of the trials from either the Music or Nature 
% experiment to train the model, and the remaining half to evaluate its performance.
%
% Finally, the script visualizes the learned weights and plots the ROC curve 
% for the trained model.
%% Add required files to path
clc
clear
close all
restoredefaultpath
addpath('Model')
addpath('visualization')
addpath('DREX-model-master')

%% Initialize
% choose the experiment type: 'nature' or 'music'
experiment_type = 'nature';


threshold = 1:-0.1:0; % threshold for computinig the ROC, the model output (alpha(t)) is compared to this threshold
fs_new = 10; % sampling rate that the features are downsampled to.

% surprisals normalization parameters, used in the "get_all_spikes" function:
normalize_threshold = 3; 
normalize_scale = 10;

%% Train and test Model on groundtruth data

% file names: name of the files generated thoughout the training process,
% please don't change
all_spikes_filename = 'Model/all_sal_new.mat';
testData_filename = 'Model/testData.mat';
timings_target_filename = 'Model/timings_target.mat';
interaction_filename = 'Model/interaction_gt';
weight_filename = 'Model/weights_gt.mat';
classification_interaction_filename = 'Model/interaction_DREX_music_gt';

% ----- To demonstrate how the model works, we first train on half of the trials -----
% We are using the groundtruth data for training in this example.

% load acoustic features for the trials we will be training on:
load(string(experiment_type)+'/all_raw_features_'+string(experiment_type)'+'_testData1.mat');

% load trial data (conditions, onset times, etc.):
load(string(experiment_type)+'/testData_block1.mat');


% "get_all_spikes" prepares features for D-REX, and computes surprisals for each feature.
% The resulting matrix "all_spikes" is the matrix that contains the normalized per-feature surprisals.
all_spikes = get_all_spikes(SoundInfoList,all_raw_features,fs_new,normalize_threshold,normalize_scale,experiment_type);
save(all_spikes_filename,'all_spikes')

% train the model and get optimized feature weights, "W"
W = train_model(SoundInfoList, testData_filename, timings_target_filename, interaction_filename, weight_filename, fs_new, [], ...
    [], experiment_type, 'gt');

% ----- testing the resulting weights using the remaining half of the trial -----



% load test data and trial info
load(string(experiment_type)+'/all_raw_features_'+string(experiment_type)+'_testData2.mat');
load(string(experiment_type)+'/testData_block2.mat');


% compute normalized surprisals for the test trials
all_spikes = get_all_spikes(SoundInfoList,all_raw_features,fs_new,normalize_threshold,normalize_scale,experiment_type);
save(all_spikes_filename,'all_spikes')

% compute boosted surprisals using the learned weights:
st_interaction_classification(convertStringsToChars(classification_interaction_filename),W,fs_new,experiment_type);
interaction1 = load(classification_interaction_filename);
all_sal = [interaction1.all_sal];

% generate testData: groundtruth data showing whether a trial is a control or if it contains a deviant. 
testData = [];
for trial_idx=1:size(all_sal,1)
    if SoundInfoList{trial_idx}.is_ctrl
        testData{trial_idx,4} = 1;
    else
        testData{trial_idx,4} = 0;
    end
end

% Compare the model's output (alpha(t)) against testData, and compute True Positive Rate, False Positive Rate and Area Under the ROC:
[TP,TN,FP,FN] = hit_miss_classification(testData, all_sal, threshold);
TP_rate_gt = TP./(TP+FN);
FP_rate_gt = FP./(FP+TN);
auc_gt = trapz(FP_rate_gt,TP_rate_gt);
acc_gt = (TP+TN)./(TP+TN+FP+FN);

%% Plot the weights and the ROC
close all
figure
draw_weights(W,'GT',experiment_type,0);
figure
plot(FP_rate_gt,TP_rate_gt,'color',[0.4660 0.6740 0.1880],'Linewidth',2)
hold on
plot(0:0.01:1,0:0.01:1,'k--')
xlim([0 1])
ylim([0 1])
xlabel('FPR')
ylabel('TPR')
title('ROC - '+string(experiment_type))
