function [save_path, keyword]= create_path(data_name,chosen_method)
% create the saving path using data name, date and chosen methods
% if there are already some saved results on the same dataset with the same
% methods on the same date, the function will increase the order number
data_path = strrep( pwd ,'\data','');
save_path = [pwd,'\result'];
date = char(datetime('today','Format','yyyy_MM_dd'));
files = dir(fullfile(save_path,'*.mat'));
exp_num = 1;
if exist(save_path,'dir') == 0
    mkdir(save_path);
end
char_chosen_method = [];
for k = 1:length(chosen_method)
    if k ~=length(chosen_method)
       char_chosen_method = [char_chosen_method,char(chosen_method(k)),'_'];
    else
       char_chosen_method = [char_chosen_method,char(chosen_method(k))];
    end
end
save_path = [save_path,'\result_','(',data_name,')_',date,'_',...
    char_chosen_method];
keyword = ['result_','(',data_name,')_',date,'_',...
    char_chosen_method]; 
for k = 1:length(files)
    if ~isempty(strfind(files(k).name,keyword))
        exp_num = exp_num + 1;
    end
end
save_path = [save_path,'_',num2str(exp_num)];

end

