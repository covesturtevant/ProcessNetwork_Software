function [C,nCounts]=getCountMat(tupleMat,nBinX,nBinY)

%This function takes only classified signals, where the signal has already
%been re-classified as integer positives representing the "bin" or
%vocabulary character that each portion of the signal is classified as.
%This will also require information as to the size of the vocabulary used
%to classify the variables. Note that there must be no missing data in
%tupleMat (i.e. rows with any missing data must have been already removed.)

nCounts = size(tupleMat,1);
C = accumarray(tupleMat,1,[nBinX nBinY nBinY]);

