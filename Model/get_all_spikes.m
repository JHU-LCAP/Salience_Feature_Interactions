function all_spikes = get_all_spikes(SoundInfoList,all_raw_features,fs_new,normalize_threshold,normalize_scale,exp)

if strcmp(exp,'music')
    exp_length = 5; % length of the experiment (s)
    base_pitch = 200; 
elseif strcmp(exp,'nature')
    exp_length = 6;
    base_pitch = 4900;
end

all_spikes_orig = cell(size(all_raw_features,2),4);
% setup D-REX parameters
params= [];
params.D = 1; % temporal dependency
params.distribution = 'gaussian'; % distribution type
params.hazard = 0.5; % hazard rate

for i=1:size(all_raw_features,2)
    for f=1:5

        if f==1 % compute surprisals for the envelope of the signal

            x = resample(all_raw_features{i}{f},fs_new,16000); % downsample feature to 10Hz
            
            params.prior = estimate_suffstat(x(1:ceil(fs_new*exp_length/4))',params); % compute the prior from the first 1/4 of the signal
            out = run_DREX_model(x',params);
            all_spikes_orig{i,4} = [all_spikes_orig{i,4} out.surprisal];

        elseif f==2 % compute surprisals for the pitch estimates

            cleaned  = clean_pitch(all_raw_features{i}{f}); % clean the pitch estimates, remove sudden large drops or rises, etc.
            x = resample(cleaned,fs_new,125); % downsample to 10Hz
            x(1:4)=x(5); % clean a bit more to remove the edge fluctuations
            x(end-1:end) = x(end-2);
            x = 2.^((2.^x)./base_pitch); % take it to exponantial scale, normalizing by the average pitch (the original pitch estimate was in log scale hence we used 2.^ twice)
           
            params.prior = estimate_suffstat(x(1:ceil(fs_new*exp_length/4)),params); % compute prior distribution
            out = run_DREX_model(x,params);
            out.surprisal(end) = out.surprisal(end-1);
            all_spikes_orig{i,4} = [all_spikes_orig{i,4} out.surprisal];

        elseif f==3 % compute surprisals for the low- and high-frequency channels of the auditory spectrogram

            x = resample(mean(all_raw_features{i}{f}(1:64,:),1),fs_new,125);
            params.prior = estimate_suffstat(x(1:ceil(fs_new*exp_length/4))',params);
            out = run_DREX_model(x',params);
            all_spikes_orig{i,4} = [all_spikes_orig{i,4} out.surprisal];

            x = resample(mean(all_raw_features{i}{f}(65:end,:),1),fs_new,125);
            params.prior = estimate_suffstat(x(1:ceil(fs_new*exp_length/4))',params);
            out = run_DREX_model(x',params);
            all_spikes_orig{i,4} = [all_spikes_orig{i,4} out.surprisal];

        elseif f>3 % compute surprisals for scale and rate
            x = resample(mean(all_raw_features{i}{f},1),fs_new,125);
            params.prior = estimate_suffstat(x(1:ceil(fs_new*exp_length/4))',params);
            out = run_DREX_model(x',params);
            all_spikes_orig{i,4} = [all_spikes_orig{i,4} out.surprisal];
        end
    end

    for f=1:6 % remove the fluctuations in the beginning and end of the computed surprisals (high values at the edges)
        all_spikes_orig{i,4}(end-1:end,f) = all_spikes_orig{i,4}(end-2,f);
        all_spikes_orig{i,4}(1,f) = all_spikes_orig{i,4}(2,f);

    end
end

all_spikes = normalize_spikes(SoundInfoList,all_spikes_orig,normalize_scale,normalize_threshold); % normalize to zero out values below the 70th percentile
end 

function cleaned  = clean_pitch(x) % function to clean the pitch using its first derivative
x(1:14)=x(15);
x(end-15:end)=x(end-15);
m = find(diff(x)>=mean(diff(x))+2*std(diff(x)));
n = find(diff(x)<=mean(diff(x))-2*std(diff(x)));
keepIndex = true(size(m));
for j = 1:length(m)-1
    % Check the difference with the next element
    diffNext = abs(m(j) - m(j + 1));

    % If the difference is less than 15, mark one of them for removal
    if diffNext<15
        keepIndex(j) = false;  % Remove the next element
    end
end
m = m(keepIndex);

keepIndex = true(size(n));
for j = 1:length(n)-1
    % Check the difference with the next element
    diffNext = abs(n(j) - n(j + 1));

    % If the difference is less than 15, mark one of them for removal
    if diffNext<15
        keepIndex(j+1) = false;  % Remove the next element
    end
end
n = n(keepIndex);
for j=1:length(n)
    [~,closest] = min(abs(m-n(j)));
    x(n(j):m(closest))=x(n(j)-1);
end
cleaned = x;
end


function all_spikes = normalize_spikes(SoundInfoList,all_spikes_orig,normalize_scale,normalize_threshold)
% first normalize between 0 and 1
maxx = [];
minn = [];
for i=1:size(SoundInfoList,1)
    maxx(i,:) = max(all_spikes_orig{i,4},[],1);
    minn(i,:) = min(all_spikes_orig{i,4},[],1);
end
max_features = max(maxx,[],1);
min_features = min(minn,[],1);
for f=1:6
    features = [];
    for i=1:size(SoundInfoList,1)
        features = [features; all_spikes_orig{i,4}(:,f)];
    end

    for i=1:size(SoundInfoList,1) % zero out the 70th percentile
        all_spikes{i,4}(:,f) = (all_spikes_orig{i,4}(:,f)-min_features(f))./(max_features(f)-min_features(f))*normalize_scale;
        all_spikes{i,4}(all_spikes{i,4}(:,f)<normalize_threshold,f) = 0;
    end
end
end
