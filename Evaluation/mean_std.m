function [m,s] = mean_std(x)
%MEAN_STD 
m = mean(x,'all');
s = std(x,0,'all');
end

