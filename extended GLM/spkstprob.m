
function [ignore_spkst,SortingProbIndex] = spkstprob(corr)
n = size(corr,2);
SortingProbIndex = nan(n,n);
ignore_spkst = eye(n,n);
for pre = 1:n-1
    for post = pre+1:n
        if isempty(corr{pre,post}) 
            continue
        end       
        mu(1) = mean(corr{pre,post}(50:52)); %center 3 bins 
        mu(2) = mean(corr{pre,post}(47:49)); %left 3 bins
        mu(3) = mean(corr{pre,post}(53:55)); %right 3 bins
        SortingProbIndex(pre,post) = nanmin( abs(mu(1)-mu(2))/mu(2), abs(mu(1)-mu(3))/mu(3) );
        SortingProbIndex(post,pre) = SortingProbIndex(pre,post);
    end
end
ignore_spkst(SortingProbIndex >= .5) = 1;
ignore_spkst(isnan(SortingProbIndex)) = 1;
