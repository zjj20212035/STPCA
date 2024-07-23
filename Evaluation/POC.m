function [POC] = POC(result, DiscFea)
%POC: Proportion of correctly selected features
L = length(result.OUTPUT);
flag = 0;
for k = 1:L
    A = result.OUTPUT{k}.FS;
    for h = 1:length(DiscFea)
        a = DiscFea(h);
        if ~isempty(find(A==a,1))
            flag = flag + 1;
        end
    end
end

POC = flag/(L*length(A));

end

