function [mean_acc,mean_nmi,std_acc,std_nmi, mean_ss] = Evaluation(data,label,id,NumFS,class_num,N,setting)
%calculate accuracy and normalized mutual information
ACC = zeros(1,N); 
NMI = zeros(1,N);
SS = zeros(1,N);    %silhouette_score
n = length(label);
data_tensor_type = setting.data_tensor_type;
% For different type of tensor data, the way of feature selection is
% different. For third-order tensor with slice-wise variables, elements in 
% the same row from one frontal slice belong to the same feature. 

    if data_tensor_type == "none" || data_tensor_type == "tube-wise" || data_tensor_type == "element-wise"
        if  ~isfield(setting,'selected_or_discarded')
            X_r = data(id(1:NumFS),:);
        elseif setting.selected_or_discarded == "selected" 
            X_r = data(id(1:NumFS),:);
        elseif setting.selected_or_discarded == "discarded" 
            X_r = data(id(NumFS+1:end),:);
        else
            error('uncorrect value of selected_or_discarded')
        end
        if ~isreal(data)
            X_r = [abs(X_r);angle(X_r)];
        end
        parfor k = 1:N
            idx = kmeans(X_r',class_num);
            final_idx = BestMapping(label,idx );
            NMI(k) = nmi(label, final_idx);
            ACC(k) = 1 - length(find(final_idx -label)) / n;
            SS(k) = mean(silhouette(X_r',label));
        end
    elseif data_tensor_type == "slice-wise"
        if  ~isfield(setting,'selected_or_discarded')
            X_r = data(id(1:NumFS),:,:);
        elseif setting.selected_or_discarded == "selected" 
            X_r = data(id(1:NumFS),:,:);
        elseif setting.selected_or_discarded == "discarded" 
            X_r = data(id(NumFS+1:end),:,:);
        else
            error('uncorrect value of selected_or_discarded')
        end
        if ~isreal(data)
            X_r = cat(1,abs(X_r),angle(X_r));
        end
        X_r = reshape(X_r,[size(X_r,1)*size(X_r,2),size(X_r,3)]);
        parfor k = 1:N
            idx = kmeans(X_r',class_num);
            final_idx = BestMapping(label,idx );
            NMI(k) = nmi(label, final_idx);
            ACC(k) = 1 - length(find(final_idx -label)) / n;
            SS(k) = mean(silhouette(X_r',label));
        end
    end

mean_acc = mean(ACC);
std_acc = std(ACC);
mean_nmi = mean(NMI);
std_nmi = std(NMI);
mean_ss = mean(SS);
end





