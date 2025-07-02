function [output] = All_Features(X,para,para_special)
%All_Features 
%No feature selection
X = X;
para_special = para_special;
nSamp = para.nSamp;
nFea = para.nFea;
output.id = 1:nSamp;
output.score = 1:nFea;
end

