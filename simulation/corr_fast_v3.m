function [d,deltaT] = corr_fast_v3(st1,st2,Ta,Tb,bin) %#codegen

% Tlist_pre : presynaptic spiking times
% Tlist_post : postsynaptic spiking times
% Ta : starting point of the histogram before the presynaptic onset
% Tb : ending interval
% bin : number of bins
% plot : 0/1 plot/don't plot

N1 = length(st1);
N2 = length(st2);

if N1 == 0 || N2 == 0
    d = [];
    deltaT = [];
    return;
end

% rough estimate of # of time difference required (assuming independence)
maxTTT = (Tb-Ta);
eN = ceil((max(N1, N2))^2 * maxTTT * 2 / min(st1(end), st2(end)));
deltaT = zeros(10 * eN, 3);

% Compute all the time differences
lastStartIdx = 1;
k = 1;
for n = 1:N1
    %     incIdx = find((st2(lastStartIdx:N2) - st1(n) >= Ta), 1, 'first');
    %     if ~isempty(incIdx)
    %         lastStartIdx = lastStartIdx + incIdx(1) - 1;
    % disp(lastStartIdx);
    incIdx = 0;
    for m = lastStartIdx:N2
        timeDiff = st2(m) - st1(n);
        if timeDiff >= Ta
            if incIdx==0
                incIdx = m;
            end
            if timeDiff <= Tb
                deltaT(k,:) = [timeDiff n m];
                k = k + 1;
            else % this is the ending point
                break;
            end
        end
    end
    if incIdx>0
        lastStartIdx = incIdx;
    end
end
deltaT = deltaT(1:(k-1),:);
edges = linspace(Ta,Tb,bin);
d = histc(deltaT(:,1),edges);