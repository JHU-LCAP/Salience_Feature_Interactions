% X is nxp matrix
% unknown c is px1
function f = detection_optimize_objective(c,X,Y)

f = 0;

% the log-likelihood loss funuction:
for i=1:size(X,1)
    f = f + log2((Y(i) + (-1).^Y(i)*(2).^(1-Y(i)).*exp(-X(i,:) * c))/(1+exp(-X(i,:) * c)));
end

f = -f;
