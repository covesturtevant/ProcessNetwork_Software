function [E] = entropyFunction(classifiedData,targetVars,lagVect,nBinMat,NoDataCode,parallelWorkers)
% Added output of linear redundancy, L
% Added selection of 'target' variables. 

[nData,nSignals] = size(classifiedData);
[nLags,~] = size(lagVect); %nLags includes the zero lag which is first
lagRange = [min(lagVect) max(lagVect)]; % 0 must be included somewhere in lagVect
nYw = 1; % the number of data points signifiying the previous history of Y. Hard-coded as 1 point previous for now, but code structured to make this an option in the future
nTuples = nData+lagRange(1)-max([lagRange(2) nYw])-1;

% INITIALIZE PARALLEL OUTPUTS OF THE SHANNON BIT FUNCTION
HXt                     = NaN(nSignals,nSignals,nLags);
HYw                     = NaN(nSignals,nSignals,nLags);
HYf                     = NaN(nSignals,nSignals,nLags);
HXtYw                   = NaN(nSignals,nSignals,nLags);
HXtYf                   = NaN(nSignals,nSignals,nLags);
HYwYf                   = NaN(nSignals,nSignals,nLags);
HXtYwYf                 = NaN(nSignals,nSignals,nLags);
I                       = NaN(nSignals,nSignals,nLags);
T                       = NaN(nSignals,nSignals,nLags);
L                       = NaN(nSignals,nSignals,nLags);
nCounts                 = NaN(nSignals,nSignals,nLags);

% PARALLELIZE ON MULTIPLE TIME LAGS
parfor (i = 1:nLags,parallelWorkers) 
                
    [   HXt(:,:,i), ...
        HYw(:,:,i), ...
        HYf(:,:,i), ...
        HXtYw(:,:,i), ...
        HXtYf(:,:,i), ...
        HYwYf(:,:,i), ...
        HXtYwYf(:,:,i), ...
        I(:,:,i), ...
        T(:,:,i), ...
        L(:,:,i), ...
        nCounts(:,:,i) ] ...
        = ShannonBitsWrapper...
        (classifiedData, ...
        targetVars, ...
        lagVect(i), ...
        nTuples, ...
        nBinMat, ...
        lagRange,...
        nYw, ...
        NoDataCode);
            
end

% ASSIGN ALL OUTPUT VARIABLES TO THE TRANSFER STRUCTURE "E"
E.HXt=HXt;
E.HYw=HYw;
E.HYf=HYf;
E.HXtYw=HXtYw;
E.HXtYf=HXtYf;
E.HYwYf=HYwYf;
E.HXtYwYf=HXtYwYf;
E.I=I;
E.T=T;
E.L=L;
E.nCounts=nCounts;



