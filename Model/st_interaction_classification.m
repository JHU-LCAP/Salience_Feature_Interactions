function st_interaction_classification(savename,W,dumb_fs,experiment)

% This function uses the optimized weights to compute the boosted surprisals (variable called "all_sal").

    % ---------------------------------------------------------------------
    % inits ---------------------------------------------------------------
    % ---------------------------------------------------------------------

    global all_spikes;
    load('all_sal_new.mat')

 
    if isempty(all_spikes)
        all_spikes = importdata('all_sal_new.mat');
    end

    if strcmp(experiment,'music')
        exp_length = 5;
    elseif strcmp(experiment,'nature')
        exp_length = 6;
    end
    range = 1;
    
    % ---------------------------------------------------------------------
    % optimization --------------------------------------------------------
    % ---------------------------------------------------------------------
    
    M_all = cell(size(all_spikes,1),1);
    for trial_idx=1:size(all_spikes,1)
        spikes = all_spikes{trial_idx,4};
        if isempty(spikes)
            continue;
        end
        M_all{trial_idx} = mat_boundary_adjust(get_M);

    end

    
    % ---------------------------------------------------------------------
    % interaction ---------------------------------------------------------
    % ---------------------------------------------------------------------
    y = [];
    M = [];
    for trial_idx=1:size(all_spikes,1)
        spikes = all_spikes{trial_idx,4};
        if isempty(spikes)
            continue;
        end
        get_interaction;
        all_spikes{trial_idx,1} = y;
    end
    all_sal = all_spikes;
    save([savename '.mat'],'all_sal');
    
    % ---------------------------------------------------------------------
    % nest functions ------------------------------------------------------
    % ---------------------------------------------------------------------
   

    function M = get_M
        M = zeros(size(spikes,1),size(spikes,2)^2);
        for t=1:size(spikes,1)
            M_sub = zeros(size(spikes,2));
            for i=1:size(spikes,2)
                for j=1:size(spikes,2)
                    if i==j
                        M_sub(i,j) = spikes(t,i);
                    else
                        % the rows are first feature
                        M_sub(i,j) = spikes(t,i) * max(spikes(max(t-range,1):min(t+range,size(spikes,1)),j));
                    end
                end
            end
            % 1-1 1-2 (pitch effect on intensity) 1-3.... 2-1 (intensity effect on pitch) 2-2 2-3 ....
            M(t,:) = reshape(M_sub',1,[]);
        end
    end

    function get_interaction
        M = get_M;
        y = zeros(size(M,1),6);
        for feat_idx=1:length(W)
            tt = M(:,(feat_idx-1)*6+(1:6)) * W(feat_idx,:)';
            y(:,feat_idx) =  2./(1+exp(-tt))-1;
        end
    end
    
    function m = mat_boundary_adjust(m)
        m = m(1:exp_length*dumb_fs,:);
    end

end


