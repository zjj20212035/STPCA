function [POC,result] = POC(result, DiscFea)
%POC: Proportion of correctly selected features
L = length(result.OUTPUT);
total = 0;
maximal_POC = 0;
maximal_POC_FS = 0;
for k = 1:L
    flag = 0;
    A = result.OUTPUT{k}.FS;
    for h = 1:length(DiscFea)
        a = DiscFea(h);
        if ~isempty(find(A==a,1))
            flag = flag + 1;
        end
    end
    total = total+flag/length(A);
    if (flag/length(A)) > maximal_POC
       maximal_POC = flag/length(A);
       maximal_POC_FS = result.OUTPUT{k}.FS;
    end
end
result.maximal_POC = maximal_POC;
result.maximal_POC_FS = maximal_POC_FS;

POC = total/L;

end

