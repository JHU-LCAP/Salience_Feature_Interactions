function M_all = st_interaction(savename,dumb_fs,experiment)

    % ---------------------------------------------------------------------
    % inits ---------------------------------------------------------------
    % ---------------------------------------------------------------------

    global testData;
    global timings_target;
    global all_spikes;
    load('testData.mat')
    load('all_sal_new.mat')
    load('timings_target.mat')
    if isempty(testData)
        disp('testData not loaded - go run et_anal first!');
        return;
    end
    if isempty(all_spikes)
        all_spikes = importdata('all_sal_new.mat');
    end

    range = 1;

    if strcmp(experiment,'music')
        target_length = 1.2;
        exp_length = 5;
    elseif strcmp(experiment,'nature')
        target_length = 1.4;
        exp_length = 6;
    end

    % ---------------------------------------------------------------------
    % optimization --------------------------------------------------------
    % ---------------------------------------------------------------------
    
    M_all = cell(size(all_spikes,1),1);
    y_all = cell(size(all_spikes,1),1);
    for trial_idx=1:size(all_spikes,1)
        spikes = all_spikes{trial_idx,4};
        if isempty(spikes)
            continue;
        end
        y_all{trial_idx} = mat_boundary_adjust(get_y);
        M_all{trial_idx} = mat_boundary_adjust(get_M);
    end
    W = get_W(cell2mat(M_all),cell2mat(y_all));
    save Model/new_W.mat W
    
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
        all_spikes{trial_idx} = y;
    end
    all_sal = all_spikes;
    save([savename '.mat'],'all_sal');
    
    % ---------------------------------------------------------------------
    % nest functions ------------------------------------------------------
    % ---------------------------------------------------------------------
    
    function y = get_y
        y = zeros(size(spikes));
        is_control = testData{trial_idx,4};
        w_target = 0.02;

        if ~is_control 
            tt = timings_target(trial_idx);
            target_t = round((tt-w_target)*dumb_fs):round((tt+target_length+w_target*2)*dumb_fs);
            target_t = target_t(target_t>0);
            condition = 1:6;
            y(target_t,condition) = ones(length(target_t),length(condition));
        end
    end

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


function W = get_W(M,Y)

    W = zeros(size(Y,2));
    % 1-2 (pitch effect on intensity) 2-1 (intensity effect on pitch) 

    W_init = [
        1 1 0 0 0 0
        1 1 1 1 0 1
        0 0 1 1 1 1
        0 0 1 1 1 1
        0 0 1 1 1 0
        0 1 1 1 0 1
    ];
    
    for feat_idx=1:6
        disp(['Running optimization for feature ' num2str(feat_idx)]);
%         tic
        
        cols = (feat_idx-1)*6 + (1:6);
        rows = find(sum(M(:,cols),2)>0);
        
        M_feat = M(rows,cols);
        Y_feat = Y(rows,feat_idx);
        
        w_feat = fmincon(@(c) detection_optimize_objective(c,M_feat,Y_feat),W_init(feat_idx,:)',...
            -eye(6),zeros(6,1),[],[],[],[],[],optimset('Algorithm','interior-point','Display','off'));
        W(feat_idx,:) = w_feat';
%         toc
    end
    
end
























