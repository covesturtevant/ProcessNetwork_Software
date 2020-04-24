function [Surrogates] = waveletVarianceMatchSurrogates(Data,Surrogates)
% Create the surrogate data for statistical testing
% 
% ----------- Inputs -----------
% Data nSamp X nVar
% Surrogates nSamp X nVar X nSur
%
% ---------- Outputs -----------
% Surrogates = the surrogate data, with the same variance at each wavelet
% timescale as those of the data
%
% ------------------------------
warning('off')

[nSamp,nVar,nSur] = size(Surrogates);
idx = 1:nSamp;

% If gaps in Data are present, fill them with linear interpolation
colGf=find(sum(isnan(Data),1) > 0);
for iVar=colGf
   
   inn=find(~isnan(Data(:,iVar)));
   Data(:,iVar) = interp1(inn,Data(inn,iVar),idx,'linear',...
       nanmean(Data(inn,iVar)));
    
end

% If gaps in Surrogates are present, fill them with linear interpolation
for iRep = 1:nSur
    colGf=find(sum(isnan(Surrogates(:,:,iRep)),1) > 0);
    for iVar=colGf

       inn=find(~isnan(Surrogates(:,iVar,iRep)));
       Surrogates(:,iVar,iRep) = interp1(inn,Surrogates(inn,iVar,iRep),...
           idx,'linear',nanmean(Surrogates(inn,iVar,iRep)));

    end
end

% How many detail scales?
detailScales = 1:floor(log2(nSamp));

% Wavelet decompose Data 
logwrite('Wavelet decomposing the data (for variance matching the surrogates)...',1);
waveData = NaN(nSamp,nVar,length(detailScales)+1);
for di = detailScales
    waveData(:,:,di) = waveletTransform(Data,di,'la8',1,1);
end
waveData(:,:,di+1) = waveletTransform(Data,di,'la8',2,1); % Approximation

% Get variances across scales for the data
stdData=std(waveData,0,1);

% Wavelet decompose the surrogates
logwrite('Wavelet decomposing & variance matching the surrogates...',1);
for iRep = 1:nSur
    logwrite(['     Surrogate ' num2str(iRep)],1);

    waveSurr = NaN(nSamp,nVar);
    for di = detailScales
        waveSurr(:,:,di) = waveletTransform(Surrogates(:,:,iRep),di,'la8',1,1);
        
        % Scale the variance
        stdSurr=std(waveSurr(:,:,di),0,1);
        waveSurr(:,:,di)=waveSurr(:,:,di).*repmat(stdData(:,:,di)./stdSurr,nSamp,1,1);

    end
    
    % Approximation
    waveSurr(:,:,di+1) = waveletTransform(Surrogates(:,:,iRep),di,'la8',2,1);
        
    % Scale the variance
    stdSurr=std(waveSurr(:,:,di+1),0,1);
    waveSurr(:,:,di+1)=waveSurr(:,:,di+1).*repmat(stdData(:,:,di+1)./stdSurr,nSamp,1,1);

    % Sum the detail and approximation to create the wavelet-variance-matched surrogate
    Surrogates(:,:,iRep) = sum(waveSurr,3);
end


warning('on')