% This function is used to find optimized feature interaction weights. 
% You need to specify whether you are training on groundtruth, behavioral
% data, or by random shuffling. The function computes the training deviant onset times
% based on the data type the model is trained on.
% This function does not include PDR training.
function W = train_model(SoundInfoList, testData_filename, timings_target_filename, interaction_filename, weight_filename, dumb_fs, pupil_info_keep, ...
    all_spikes, experiment, type_exp)

% specify experiment length, used when training on PDR:
if strcmp(experiment,'music')
    exp_length = 5000; % in samples, 1kHz sampling frequency for the PDR data
elseif strcmp(experiment,'nature')
    exp_length = 6000;
end

% Initialize
testData = [];
timings_target = [];
num_trials = size(SoundInfoList,1);

% --------- Training on Ground Truth (GT) ---------
if strcmp(type_exp,'gt')
    for i=1:num_trials
        if SoundInfoList{i}.is_ctrl % no change reported (control trial)
            testData{i,4} = 1; % Mark as control
            timings_target(i) = 0; % no deviant onset time
        else
            testData{i,4} = 0; % mark as non-control, containing deviant
            timings_target(i) = SoundInfoList{i}.fg_onset_time; % Use known foreground deviant onset time
        end
    end

% --------- Training on Behavioral Data ---------   
elseif strcmp(type_exp,'beh')
    for i=1:num_trials
        behavior_idx = find(pupil_info_keep(:,13)==i); % Find trial index in pupil data

        % Response code 3 or 4 = "no change" or "miss"
        if pupil_info_keep(behavior_idx,9)==3 || pupil_info_keep(behavior_idx,9)==4 % no change reported  % pupil_info(:,9): subject response
            testData{i,4} = 1;
            timings_target(i) = 0;
        else
            testData{i,4} = 0;
            timings_target(i) = SoundInfoList{i}.fg_onset_time;
            if pupil_info_keep(behavior_idx,9)==1 % true positive: use known groundtruth onset time
                timings_target(i) = SoundInfoList{i}.fg_onset_time;
            elseif pupil_info_keep(behavior_idx,9)==2 % false positive: Estimate onset based on surprisal data 
                y = sum(all_spikes{i,4},2)/6;
                [m,n] = find(y==max(y,[],'all')); % find maximum sum of surprisals
                timings_target(i) = round((m(1))/dumb_fs,1); % Convert to seconds
            end
        end

    end

% --------- Training on Randomized Labels (for comparison with the behavioral model)---------
elseif strcmp(type_exp,'rand')
    for i=1:num_trials
        % use behavioral data to setup deviant onset times
        behavior_idx = find(pupil_info_keep(:,13)==i);
        if pupil_info_keep(behavior_idx,9)==3 || pupil_info_keep(behavior_idx,9)==4 % no change reported
            testData{i,4} = 1; % pupil_info(:,9): subject response
            timings_target(i) = 0;
        else
            testData{i,4} = 0;
            timings_target(i) = SoundInfoList{i}.fg_onset_time;
            if pupil_info_keep(behavior_idx,9)==1 % true positive
                timings_target(i) = SoundInfoList{i}.fg_onset_time;
            elseif pupil_info_keep(behavior_idx,9)==2 % false positive
                y = sum(all_spikes{i,4},2)/6;
                [m,n] = find(y==max(y,[],'all'));
                timings_target(i) = round((m(1))/dumb_fs,1);
            end
        end

    end

    % Shuffle trial labels and timings randomly
    shuffled_idx = randperm(size(testData,1));
    testData_hold = [];
    for i=1:size(testData,1)
        testData_hold{i,4} = testData{shuffled_idx(i),4};
    end
    testData = testData_hold;
    timings_target = timings_target(shuffled_idx);
end
save(testData_filename,'testData')
save(timings_target_filename, 'timings_target')
str_name = convertStringsToChars(interaction_filename);

% Run interaction analysis script to get interaction weights
st_interaction(str_name,dumb_fs,experiment);
new_W = load('new_W.mat');
W = new_W.W;
save(weight_filename, 'W')
end

