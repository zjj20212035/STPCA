function [output,acc,nmi,ss] = ClusterExp(struct,result)
% struct: a struct that contains the dataset and all parameters.
% struct.X: data matrix of nFea*nSamp (nFea:#Feautre;nSamp:#Sample)
%    or data tensor with the last dimension being 'nSamp' 
% struct.label: just labels
% struct.method: a handle of the function of certain method. @method
% struct.Ten_based: whether the method is tensor-based or not (0:No; 1:Yes)
% struct.tensor_type: the type of data tensor: "tube-wise",or "slice-wise"
% struct.tensor_size: the size of a tensor sample
% struct.Dname: name of the dataset
% struct.Mname: name of the method
% struct.para: parameters needed for grid search in @method (a struct with size [1,#parameters])
% struct.para_special: parameters that are not for grid search in @method
% struct.NumFS: number of selected features (can be a scalar or a vector)
% struct.tensor_size: size of a single input tensor (e.p. image size or video size)
% struct.figpara: parameters needed for visualization of feature selection
% struct.N: repeated times of clustering
% struct.save_path: path for saving
data = struct.data;
label = struct.label;
if ~isfield(struct,'tensor_type')
    data_tensor_type = "none";
else
    data_tensor_type = struct.tensor_type;
end
Dname = struct.Dname;
Fname = struct.Fname;
para = struct.para;
NumFS = struct.NumFS;
if isfield(struct,'tensor_type')
    tensor_size = struct.tensor_size;
end
N = struct.N;
[para_grid, para_grid_struct, size_grid, para_fieldNames]= GenGrid(para);
save_path = [struct.save_path,'_temp'];


%% handle the dataset
size_X = size(data);
if length(size_X) == 2
   data_vec = data;
else  
   nSamp = size_X(end);
   nFea = prod(size_X(1:end-1));
   data_vec = reshape(data,[nFea,nSamp]);
end   
NumClass= length(unique(label));  

%% create a series of tensors or cells to save the results
L_NumFS = length(NumFS);
ACC = tenzeros([L_NumFS,size_grid]); % to save accuracy
STD_ACC = tenzeros([L_NumFS,size_grid]); % to save standard deviation for acc
NMI = tenzeros([L_NumFS,size_grid]); % to save normalized mutual information
STD_NMI = tenzeros([L_NumFS,size_grid]);% to save standard deviation for nmi
SS = tenzeros([L_NumFS,size_grid]); %to save silhouette_score
optimal_para_acc = cell(1,L_NumFS); % to save the parameters corresponding to the max acc
optimal_para_nmi = cell(1,L_NumFS); % to save the parameters corresponding to the max nmi
optimal_para_acc_std = cell(1,L_NumFS); % to save the parameters corresponding to std of the max acc
optimal_para_nmi_std = cell(1,L_NumFS); % to save the parameters corresponding to std of the max nmi
optimal_para_ss = cell(1,L_NumFS); 

%% Clustering experiment using k-means
%%%% unfolding the metric tensor for the convenience of grid search
ACC = tenmat(ACC,1);
STD_ACC = tenmat(STD_ACC,1);
NMI = tenmat(NMI,1);
STD_NMI = tenmat(STD_NMI,1);
SS = tenmat(SS,1);
for k = 1:size(para_grid,1)
    score = result.OUTPUT{k}.score;
    id = result.OUTPUT{k}.id;
    if data_tensor_type == "none" || data_tensor_type == "tube-wise" || data_tensor_type == "element-wise"
        setting.data_tensor_type = data_tensor_type;
        if isfield(setting,'selected_or_discarded')
            setting.selected_or_discarded = selected_or_discarded;
        end
        
        parfor p = 1:length(NumFS)
            [ACC(p,k),NMI(p,k),STD_ACC(p,k),STD_NMI(p,k),SS(p,k)] = Evaluation(data_vec,label,id,NumFS(p),NumClass,N,setting);
        end
        
        
    elseif data_tensor_type == "slice-wise"
        setting.data_tensor_type = data_tensor_type;
        tensor_score = reshape(score,tensor_size);
        tensor_score = sum(tensor_score,2);
        [~,id] = sort(tensor_score,'descend');
        result.OUTPUT{k}.id = id;
        if isfield(setting,'selected_or_discarded')
            setting.selected_or_discarded = selected_or_discarded;
        end
        parfor p = 1:length(NumFS)
            [ACC(p,k),NMI(p,k),STD_ACC(p,k),STD_NMI(p,k),SS(p,k)] = Evaluation(data,label,id,NumFS(p),NumClass,N,setting);
        end
    end
    result.OUTPUT{k}.acc = (ACC(:,k))';
    result.OUTPUT{k}.nmi = (NMI(:,k))';
    result.OUTPUT{k}.std_acc = (STD_ACC(:,k))';
    result.OUTPUT{k}.std_nmi = (STD_NMI(:,k))';
    result.OUTPUT{k}.ss = (SS(:,k))';

    %% print the process
    disp([Dname,', ',Fname,', ','Clustering.',num2str(k)])
    save(save_path)
end

%%  sort results
%%%% find the max acc, the max nmi and their corresponding standard deviation
if  size(double(ACC),2) > 1
    [acc, optimal_para_position_acc] = max(double(ACC'));
    [nmi, optimal_para_position_nmi] = max(double(NMI'));
    [ss, optimal_para_position_ss] = max(double(SS'));
else
    acc = double(ACC');
    optimal_para_position_acc = ones(1,L_NumFS);
    nmi = double(NMI');
    optimal_para_position_nmi = ones(1,L_NumFS);
    ss = double(SS');
    optimal_para_position_ss = ones(1,L_NumFS);    
end

% for display
acc_disp = acc*100; nmi_disp = nmi*100;
acc_percent = arrayfun(@(x) sprintf('%.2f%%', x), acc_disp, 'UniformOutput', false);
nmi_percent = arrayfun(@(x) sprintf('%.2f%%', x), nmi_disp, 'UniformOutput', false);
acc_percent_disp = strjoin(acc_percent, ' ');
nmi_percent_disp = strjoin(nmi_percent, ' ');
ss_disp = sprintf('%.4f  ', ss); 
ss_disp = ['[', strtrim(ss_disp), ']'];
[best_acc,idx_best_acc] = max(acc_disp, [], 'all');
[best_nmi,~] = max(nmi_disp, [], 'all');
best_acc_disp = sprintf('%.2f%%', best_acc);
best_nmi_disp = sprintf('%.2f%%', best_nmi);

for k = 1:L_NumFS
    for p = 1:length(para_fieldNames)
        fieldName = para_fieldNames{p};
        optimal_para_acc{k}.(fieldName) = para_grid(optimal_para_position_acc(k),p);
        optimal_para_nmi{k}.(fieldName) = para_grid(optimal_para_position_nmi(k),p);
        optimal_para_acc_std{k}.(fieldName) = para_grid(optimal_para_position_acc(k),p);
        optimal_para_nmi_std{k}.(fieldName) = para_grid(optimal_para_position_nmi(k),p);
        optimal_para_ss{k}.(fieldName) = para_grid(optimal_para_position_ss(k),p);
    end
end

disp(['for ', Fname, ', ACC = [', acc_percent_disp, '] with [', num2str(NumFS), '] features']);
disp(['for ', Fname, ', NMI = [', nmi_percent_disp, '] with [', num2str(NumFS), '] features']);
disp(['for ', Fname, ', SS = ', ss_disp, ' with ', num2str(NumFS), '] features']);
disp(['On ',Dname, ', the best results of ', Fname, ': ACC = ', best_acc_disp, '  NMI = ', best_nmi_disp, '  SS = ', num2str(max(ss, [], 'all'))]);
%%%% fold all metric matrix into tensor for the convenience of using them again
% the first dimension is always the number of selected features
output = result;
output.ACC = tensor(ACC);
output.NMI = tensor(NMI);
output.SS = tensor(SS);
output.STD_ACC = tensor(STD_ACC);
output.STD_NMI = tensor(STD_NMI);
output.optimal_para_acc = optimal_para_acc;
output.optimal_para_nmi = optimal_para_nmi;
output.optimal_para_position_acc = optimal_para_position_acc;
output.optimal_para_position_nmi = optimal_para_position_nmi;
output.best_acc_disp = best_acc_disp;
output.best_nmi_disp = best_nmi_disp;
output.best_ss_disp = num2str(max(ss, [], 'all'));
output.time_disp = double(result.TIME(optimal_para_position_acc(idx_best_acc)));
delete([save_path,'.mat'])
end


function [para_grid, para_grid_struct, para_size, para_fieldNames]= GenGrid(para)
% to generate a parameter-grid
para_fieldNames = fieldnames(para);
NumPara = length(para_fieldNames);
para_size = zeros(1,NumPara);
strParaName = cell(1,NumPara);
for k = 1:NumPara
    fieldName = para_fieldNames{k};
    if k~=NumPara
        strParaName{k} = ['para.',fieldName,','];
    else
        strParaName{k} = ['para.',fieldName];
    end
    para_size(k) = length(para.(fieldName));
end
eval(['para_grid = (combvec(',cell2mat(strParaName),'))'';']);
length_grid = size(para_grid,1);
para_grid_struct = cell(length_grid,1);
for k = 1:length_grid
    for p = 1:NumPara
        fieldName = para_fieldNames{p};
        para_grid_struct{k}.(fieldName) = para_grid(k,p);
    end
end

end