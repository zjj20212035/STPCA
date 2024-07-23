function [POTC] = POTC(result, DiscFea)
%POTC: Proportion of the times of selecting all correct features
L = length(result.OUTPUT);
flag = 0;
for k = 1:L
    A = result.OUTPUT{k}.FS;
    A = sort(A,'ascend');
    DiscFea = sort(DiscFea,'ascend');
    if A(:) == DiscFea(:)
        flag = flag + 1;
    end
end

POTC = flag/length(result.OUTPUT);

end