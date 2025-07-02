%%%% A reproduction of Figure 7 (Visualization) in 
%%%% Orientation-Aware Sparse Tensor PCA for Efficient Unsupervised Feature
%%%% Selection%

%% preparing environment
clear
% include tool path
addpath(genpath([pwd,'\methods']))
addpath(genpath([pwd,'\Evaluation']))
addpath(genpath([pwd,'\utils']))
% start parallel computation toolbox
try 
    parpool;
catch
    disp('parpool already available')
end

%% experiment
%%%%%%%%%%%% JAFFE %%%%%%%%%%%%%%%%%
data_name = 'JAFFE';
data_path = fullfile('..', 'data', data_name);
config_Figure7;
% conduct STPCA-MP
para_special_STPCA_MP_JAFFE = para_special_STPCA_MP;
para_special_STPCA_MP_JAFFE.orientation = 1;
struct_STPCA_MP_JAFFE = struct_STPCA_MP;
struct_STPCA_MP_JAFFE.NumFS = 50:50:300;
struct_STPCA_MP_JAFFE.para_special = para_special_STPCA_MP_JAFFE;
[result_STPCA_MP_JAFFE] = AlgExecution(struct_STPCA_MP_JAFFE);
[result_STPCA_MP_JAFFE,~,nmi_STPCA_MP_JAFFE,~,...,
    ] = ClusterExp(struct_STPCA_MP_JAFFE,result_STPCA_MP_JAFFE);
ticks = [1,20:20:60];
num = 300;
[~,fea_pos_STPCA_MP_JAFFE] = max(nmi_STPCA_MP_JAFFE);
FS_STPCA_MP_JAFFE_pos = result_STPCA_MP_JAFFE.optimal_para_position_nmi(fea_pos_STPCA_MP_JAFFE);
FS_STPCA_MP_JAFFE = result_STPCA_MP_JAFFE.OUTPUT{FS_STPCA_MP_JAFFE_pos}.id;
% original image
ratio = 0.7;
id_imag = 115;
imag = data(:,:,id_imag);
figure;imshow(imag*ratio);
colormap summer
% score map and highlighted selected features of STPCA-MP
output = result_STPCA_MP_JAFFE.OUTPUT{1,FS_STPCA_MP_JAFFE_pos};
score = output.score;
score = reshape(score,tensor_size);
if max(score,[],'all')>10^-4
   score = score./max(score,[],'all');
else
   score = score*10^4;
   score = score./max(score,[],'all');
end
h=figure;
imagesc(score);
set(gca,'FontSize',25)
set(h, 'Position', [867, 337, 560, 490]);
xticks(ticks);
yticks(ticks);
colorbar
figure;
imag_STPCA_MP_Dir1 = imag*ratio;
imag_STPCA_MP_Dir1(FS_STPCA_MP_JAFFE(1:num)) = 1;
imshow(imag_STPCA_MP_Dir1);
colormap summer

%%%%%%%%%%%% PIE %%%%%%%%%%%%%%%%%
data_name = 'PIE';
data_path = fullfile('..', 'data', data_name);
config_Figure7;
% conduct STPCA-MP
para_special_STPCA_MP_PIE = para_special_STPCA_MP;
para_special_STPCA_MP_PIE.orientation = 1;
struct_STPCA_MP_PIE = struct_STPCA_MP;
struct_STPCA_MP_PIE.NumFS = 50:50:300;
struct_STPCA_MP_PIE.para_special = para_special_STPCA_MP_PIE;
[result_STPCA_MP_PIE] = AlgExecution(struct_STPCA_MP_PIE);
[result_STPCA_MP_PIE,~,nmi_STPCA_MP_PIE,~,...,
    ] = ClusterExp(struct_STPCA_MP_PIE,result_STPCA_MP_PIE);
ticks = [1,10:10:30];
num = 100;
[~,fea_pos_STPCA_MP_PIE] = max(nmi_STPCA_MP_PIE);
FS_STPCA_MP_PIE_pos = result_STPCA_MP_PIE.optimal_para_position_nmi(fea_pos_STPCA_MP_PIE);
FS_STPCA_MP_PIE = result_STPCA_MP_PIE.OUTPUT{FS_STPCA_MP_PIE_pos}.id;
% original image
ratio = 0.7;
id_imag = 50;
imag = data(:,:,id_imag);
figure;imshow(imag*ratio);
colormap summer
% score map and highlighted selected features of STPCA-MP
output = result_STPCA_MP_PIE.OUTPUT{1,FS_STPCA_MP_PIE_pos};
score = output.score;
score = reshape(score,tensor_size);
if max(score,[],'all')>10^-4
   score = score./max(score,[],'all');
else
   score = score*10^4;
   score = score./max(score,[],'all');
end
h=figure;
imagesc(score);
set(gca,'FontSize',25)
set(h, 'Position', [867, 337, 560, 490]);
xticks(ticks);
yticks(ticks);
colorbar
figure;
imag_STPCA_MP_Dir1 = imag*ratio;
imag_STPCA_MP_Dir1(FS_STPCA_MP_PIE(1:num)) = 1;
imshow(imag_STPCA_MP_Dir1);
colormap summer

%%%%%%%%%%%% BreastMNIST %%%%%%%%%%%%%%%%%
data_name = 'BreastMNIST';
data_path = fullfile('..', 'data', data_name);
config_Figure7;
% conduct STPCA-MP
para_special_STPCA_MP_BreastMNIST = para_special_STPCA_MP;
para_special_STPCA_MP_BreastMNIST.orientation = 2;
struct_STPCA_MP_BreastMNIST = struct_STPCA_MP;
struct_STPCA_MP_BreastMNIST.NumFS = 50:50:300;
struct_STPCA_MP_BreastMNIST.para_special = para_special_STPCA_MP_BreastMNIST;
[result_STPCA_MP_BreastMNIST] = AlgExecution(struct_STPCA_MP_BreastMNIST);
[result_STPCA_MP_BreastMNIST,~,nmi_STPCA_MP_BreastMNIST,~,...,
    ] = ClusterExp(struct_STPCA_MP_BreastMNIST,result_STPCA_MP_BreastMNIST);
ticks = [1,10:10:30];
num = 50;
[~,fea_pos_STPCA_MP_BreastMNIST] = max(nmi_STPCA_MP_BreastMNIST);
FS_STPCA_MP_BreastMNIST_pos = result_STPCA_MP_BreastMNIST.optimal_para_position_nmi(fea_pos_STPCA_MP_BreastMNIST);
FS_STPCA_MP_BreastMNIST = result_STPCA_MP_BreastMNIST.OUTPUT{FS_STPCA_MP_BreastMNIST_pos}.id;
% original image
ratio = 0.7;
id_imag = 79;
imag = data(:,:,id_imag);
figure;imshow(imag*ratio);
colormap summer
% score map and highlighted selected features of STPCA-MP
output = result_STPCA_MP_BreastMNIST.OUTPUT{1,FS_STPCA_MP_BreastMNIST_pos};
score = output.score;
score = reshape(score,tensor_size);
if max(score,[],'all')>10^-4
   score = score./max(score,[],'all');
else
   score = score*10^4;
   score = score./max(score,[],'all');
end
h=figure;
imagesc(score);
set(gca,'FontSize',25)
set(h, 'Position', [867, 337, 560, 490]);
xticks(ticks);
yticks(ticks);
colorbar
figure;
imag_STPCA_MP_Dir1 = imag*ratio;
imag_STPCA_MP_Dir1(FS_STPCA_MP_BreastMNIST(1:num)) = 1;
imshow(imag_STPCA_MP_Dir1);
colormap summer