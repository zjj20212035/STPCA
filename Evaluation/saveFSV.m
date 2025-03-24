function saveFSV(X,Para)
% To save the visualization results for feature selection
% X: the data matrix of size nFea*nSamp (or tensor of size m*n*nSamp)
% Para: some parameters that determine the configurations of saving
% Para.name: name of the dataset
% Para.Mname: name of the method
% Para.idGrid: index of grid search
% Para.idFea: index of selected features.
% Para.NumFS: number of the selected features
% Para.imag_size: the size of one image
% Para.type: the saving types of figures (fig or png or jpg, etc..)
%% parameter setting
Dname = Para.Dname;
Mname = Para.Mname;
idGrid = Para.idGrid;
idFea = Para.idFea;
NumFS_fig = Para.NumFS_fig;
imag_size = Para.imag_size;
type = Para.type;

%% handling X
size_X = size(X);
nSamp = size_X(end);
if length(size_X) == 3
    X = reshape(X,[prod(imag_size),size_X(2)]);
end

%% create a saving path
for k = 1:nSamp
    if exist([pwd,'\result\',Dname,'\',Mname,'\',num2str(idGrid),'\',num2str(k)]) == 0
        mkdir([pwd,'\result\',Dname,'\',Mname,'\',num2str(idGrid),'\',num2str(k)])
    end
end

%% start saving
for p = NumFS_fig
    imags = X;
    imags(idFea(1:p),:) = 1;
    parfor k = 1:nSamp
        save_path = [pwd,'\result\',Dname,'\',Mname,'\',num2str(idGrid),'\',num2str(k),...
              '\',num2str(p),'.'];
        imag = reshape(imags(:,k),imag_size);
        if type ~= "fig"
            imwrite(imag,[save_path,type]);
        else           
            figure;
            imshow(imag)
            savefig([save_path,'fig'])
            close
        end
    end
end

%% 
end

