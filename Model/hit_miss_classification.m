function [TP,TN,FP,FN] = hit_miss_classification(testData, all_sal, threshold)
TP = zeros(1,length(threshold));
TN = zeros(1,length(threshold));
FP = zeros(1,length(threshold));
FN = zeros(1,length(threshold));
for trial_idx=1:size(all_sal,1)
    for th=1:length(threshold)
        model_output = sum(all_sal{trial_idx,1}(1:end,:),2)/6; % This is alpha(t): the average over boosted surprisals
        output_spike = find(max(model_output)>threshold(th));

            if testData{trial_idx,4} && isempty(output_spike)
                TN(th) = TN(th)+1;
            end
            if testData{trial_idx,4} && ~isempty(output_spike)
                FP(th) = FP(th)+1;
            end
            if ~testData{trial_idx,4} && isempty(output_spike)
                FN(th) = FN(th)+1;
            end
            if ~testData{trial_idx,4} && ~isempty(output_spike)
                TP(th) = TP(th)+1;
            end
    end
end
end