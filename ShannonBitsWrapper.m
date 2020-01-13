function [HXt, HYw, HYf, HXtYw, HXtYf, HYwYf, HXtYwYf, I, T, L, nCounts] ...
    = ShannonBitsWrapper(classifiedData, targetVars, lag, nTuples, nBinMat, lagRange, nYw, NoDataCode)
% Added calculation and output of linear redundancy L according to Palus 2008
% Constrained looping over target Y variables only

[~,nSignals] = size(classifiedData);
classifiedData(classifiedData == NoDataCode) = NaN;

HXt             = NaN( nSignals, nSignals );
HYw             = NaN( nSignals, nSignals );
HYf             = NaN( nSignals, nSignals );
HXtYw           = NaN( nSignals, nSignals );
HXtYf           = NaN( nSignals, nSignals );
HYwYf           = NaN( nSignals, nSignals );
HXtYwYf         = NaN( nSignals, nSignals );
L               = NaN( nSignals, nSignals );
nCounts         = NaN( nSignals, nSignals );

%INITIALIZE TIME-SHIFTED COLUMNS
XtSTART = max([lagRange(2) nYw])+1-lag;
YwSTART = max([lagRange(2) nYw])+1-nYw;
YfSTART = max([lagRange(2) nYw])+1; 

for sY=targetVars
    tupleMatYw=classifiedData(YwSTART:YwSTART+nTuples-1,sY);
    tupleMatYf=classifiedData(YfSTART:YfSTART+nTuples-1,sY);
    
    for sX=1:nSignals

        %CONSTRUCT THREE-COLUMN MATRIX WITH COLUMNS TIME-SHIFTED
        %tupleMat=NaN(nTuples,3);
        tupleMat=[classifiedData(XtSTART:XtSTART+nTuples-1,sX) tupleMatYw tupleMatYf];
        %tupleMat(:,1)=classifiedData(XtSTART:XtSTART+nTuples-1,sX);        %Leading Node Xt (lag tau earlier than present)
        %tupleMat(:,2)=classifiedData(YwSTART:YwSTART+nTuples-1,sY);        %Led Node Yw (present time)
        %tupleMat(:,3)=classifiedData(YfSTART:YfSTART+nTuples-1,sY);        %Led Node Yf (one timestep in future)
        tupleMat(isnan(sum(tupleMat,2)),:) = []; % Remove an rows with missing data. No worries about destroying time relationships at this point
       
        %GET COUNTS
        [C, nCounts(sX,sY)] = getCountMat( tupleMat, nBinMat(sX), nBinMat(sY));
        
        %CHECK TO ENSURE TUPLEMAT HAS AT LEAST ONE COMPLETE ROW OF DATA
        if nCounts(sX,sY) == 0
            logwrite(['Warning: no data in tupleMat, skipping sX = ' num2str(sX) ', sY = ' num2str(sY) ', lag = ' num2str(lag)],1)
            continue
        end
 
        %CALCULATE ENTROPIES FROM TUPLEMAT
        [HXt(sX,sY),HYw(sX,sY),HYf(sX,sY),HXtYw(sX,sY),HXtYf(sX,sY),HYwYf(sX,sY),HXtYwYf(sX,sY)]=GetShannonBits( C, nCounts(sX,sY) );

        % Calculate the linear redundancy according to Palus 2008
        CM = corrcoef(tupleMat(:,[1 3])); % Correlation matrix between Xt and Yf. Requires no missing data
        CM(isnan(CM)) = 0;
        L(sX,sY) = -0.5*(sum(log2(eig(CM))));
        
    end
end

%CALCULATE ADDITIONAL ENTROPIES
I = HXt + HYf - HXtYf;
T = HXtYw + HYwYf - HYw - HXtYwYf;