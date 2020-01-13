clear

sitename = 'NZKop';
surrogateWaveDir = 'C:\Users\csturtevant\Dropbox\KnoxLab\CoveAnalysis\SensitivityTests\SampleSize\WaveletDecomposed\';
surrogateWaveDetailPrefix = 'RandomWalk100_15years_halfHourly_detail_';
surrogateWaveApproxPrefix = 'RandomWalk100_15years_halfHourly_approx_';
dataFile = 'C:\Users\csturtevant\Dropbox\KnoxLab\CoveAnalysis\testSaraData\raw\NZKop_DailyAvg_GF.mat'; % No gaps
outFile = 'C:\Users\csturtevant\Dropbox\KnoxLab\CoveAnalysis\testSaraData\raw\NZKop_DailyAvg_GF_RWvarMatchSurr.mat';

detailScales = 1:12; % 1:N
numSurr = 100; % Must match what's in the surrogate files

% Wavelet decompose the data at the scales indicated (plus the
% approximation at the final scale)
load(dataFile) % Variable Data
waveData = repmat(NaN*Data,1,1,length(detailScales)+1);
for di = detailScales
    waveData(:,:,di) = waveletTransform(Data,di,'la8',1,1);
end
waveData(:,:,di+1) = waveletTransform(Data,di,'la8',2,1);
numData = size(waveData,1);
numVar = size(waveData,2);

% Now go through the surrogate scales, matching the variance at each scale
% with that of the data
Surrogates = zeros(numData,numVar,numSurr);
for di = [detailScales detailScales(end)+1]
    
    % Open the surrogates file at this scale - variable Data
    if di <= detailScales(end)
        load([surrogateWaveDir surrogateWaveDetailPrefix num2str(di) '.mat'])
    else
        load([surrogateWaveDir surrogateWaveApproxPrefix num2str(detailScales(end)) '.mat'])
    end
    
    % Reformat the surrogates to match the data
    %numDataSurrogates = size(Data,1);
    %setSurrogatesMin = max([floor((numDataSurrogates-numData)/2) 1]);
    %setSurrogates = setSurrogatesMin:((setSurrogatesMin+numData)-1);
    %SurrVarMatch = Data(setSurrogates,:,:);
    %SurrVarMatch = repmat(reshape(SurrVarMatch,size(SurrVarMatch,1),1,size(SurrVarMatch,2)),1,size(waveData,2));
    SurrVarMatch = NaN(numData,numVar,numSurr);
    for iVar = 1:numVar
        setSurrogates = 100001+(((iVar-1)*numData+1):iVar*numData);
        SurrVarMatch(:,iVar,:) = reshape(Data(setSurrogates,:),numData,1,numSurr);
    end
    
    % Get Variances
    stdData=std(waveData(:,:,di),0,1,'omitnan');
    stdSurr=std(SurrVarMatch,0,1,'omitnan');
    
    % Variance match
    SurrVarMatch=SurrVarMatch./repmat(stdSurr,numData,1,1).*repmat(stdData,numData,1,size(SurrVarMatch,3));
    Surrogates = Surrogates + SurrVarMatch; % Add the variance-matched detail/approx at this scale
    
    clear Data
    
end

% Reload the data file
load(dataFile) % Variable Data
save(outFile,'Data','Surrogates')

