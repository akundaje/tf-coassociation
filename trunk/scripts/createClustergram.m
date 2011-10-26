signal = importdata(inputFile); % Read file
signal.data = signal.data + rand(size(signal.data))*1e-7; % Regularize data
signal.tfids = signal.textdata(1,2:end);
signal.peakids = signal.textdata(2:end,1);
signal = rmfield(signal,'textdata');
signal = dataset({signal.data,signal.tfids{:}} , 'ObsNames' , signal.peakids);
stdVals = datasetfun(@std,signal,'UniformOutput',true,'DatasetOutput',false);
signal(:,(stdVals<1e-3)) = [];

load('cmap.mat','cmap');

% create clustergram
cGram = clustergram(double(signal) , 'Standardize' , 'none' , ...
    'RowPDist' , 'correlation' , 'ColumnPDist' , 'cosine', ...
    'OptimalLeafOrder' , true , 'symmetric' , false , 'Colormap' , cmap , ...
    'RowLabels' , signal.Properties.ObsNames , 'ColumnLabels' , signal.Properties.VarNames );
addXLabel(cGram,'TFs');
[~,titleName,~] = fileparts(inputFile); 
addTitle(cGram,strrep(titleName,'_','-'));
%save([inputFile,'.data.mat']);