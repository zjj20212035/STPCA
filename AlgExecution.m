function [result] = AlgExecution(struct)
% struct: a struct that contains the dataset and all parameters.
% struct.X: data matrix of nFea*nSamp (nFea:#Feautre;nSamp:#Sample)
%    or data tensor with the last dimension being 'nSamp' 
% struct.label: just labels
% struct.method: a handle of the function of certain method. @method
% struct.Ten_based: whether the method is tensor-based or not (0:No; 1:Yes)
% struct.tensor_type: the type of data tensor: 
%        "tube-wise",or "slice-wise"
% struct.tensor_size: the size of a tensor sample
% struct.Dname: name of the dataset
% struct.Mname: name of the method
% struct.para: parameters needed for grid search in @method (a struct with size [1,#parameters])
% struct.para_special: parameters that are not for grid search in @method
% struct.NumFS: number of selected features (can be a scalar or a vector)
% struct.tensor_size: size of a single input tensor (e.p. image size or video size)
% struct.save_path: path for saving
%% unzip 'struct'
data = struct.data;
method = struct.method;
Ten_based = struct.Ten_based;
Dname = struct.Dname;
Fname = struct.Fname;
para = struct.para;
para_special = struct.para_special;
NumFS = struct.NumFS;
data_tensor_based = struct.data_tensor_based;
if data_tensor_based == 1
    tensor_type = struct.tensor_type;
end
if struct.data_tensor_based == 1
    tensor_size = struct.tensor_size;
else
    tensor_size = [size(data,1),1];
end
save_path = [struct.save_path,'_temp'];

%% handle the dataset
size_X = size(data);
if length(size_X) == 2
   nSamp = size_X(2);
   data_ten = reshape(data,[tensor_size,nSamp]);
   data_vec = data;
else  
   nSamp = size_X(end);
   nFea = prod(size_X(1:end-1));
   data_ten = data;
   data_vec = reshape(data,[nFea,nSamp]);
end   

%% create a series of tensors or cells to save the results
%%%% to generate a parameter-grid
[para_grid, para_grid_struct, size_grid, para_fieldNames]= GenGrid(para);
NumPara = length(para_fieldNames);
TIME = tenzeros([1,size_grid]); % to save training time 
OUTPUT = cell(1,prod(size_grid)); % to save the rank of feature scores

for k = 1:prod(size_grid)
    %% execute the method
    if Ten_based == 0 % for non-tensor-based methods
        tic
        output = method(data_vec,para_grid_struct{k},para_special);
        TIME(k) = toc;
        score = output.score;
        if data_tensor_based == 1
            if tensor_type == "slice-wise"
                tensor_score = reshape(score,tensor_size);
                tensor_score = sum(tensor_score,2);
                [~,id] = sort(tensor_score,'descend');
                output.id = id;
                output.FS = id(1:NumFS);
            else
                [~,id] = sort(score,'descend');
                output.id = id;
                output.FS = id(1:NumFS);
            end
        else
            [~,id] = sort(score,'descend');
            output.id = id;
            output.FS = id(1:NumFS);
        end

        output.time = TIME(k);
        OUTPUT{k} = output;
    else % for tensor-based methods
        tic
        output = method(data_ten,para_grid_struct{k},para_special);
        TIME(k) = toc;
        score = output.score;
        if data_tensor_based == 1
            if tensor_type == "slice-wise"
                tensor_score = reshape(score,tensor_size);
                tensor_score = sum(tensor_score,2);
                [~,id] = sort(tensor_score,'descend');
                output.id = id;
                output.FS = id(1:NumFS);
            else
                [~,id] = sort(score,'descend');
                output.id = id;
                output.FS = id(1:NumFS);
            end
        else
            [~,id] = sort(score,'descend');
            output.id = id;
            output.FS = id(1:NumFS);
        end
        output.time = TIME(k);
        OUTPUT{k} = output;
    end
    
    %% print the process
    disp([Dname,', ',Fname,', ','Grid searching.',num2str(k)])
    dispStr = [];
    for p = 1:NumPara
        fieldName = para_fieldNames{p};
        dispStr = [dispStr,fieldName, ' = ',num2str(para_grid_struct{k}.(fieldName)),'; '];
    end
    disp([dispStr,'training time = ',num2str(TIME(k)),'s'])
    save(save_path)

end
result.TIME = tensor(TIME);
result.OUTPUT = OUTPUT;
result.size_grid = size_grid;
result.para_grid = para_grid;
result.para_fieldNames = para_fieldNames;
save(save_path)
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

