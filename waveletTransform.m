function [waveData] = waveletTransform(rawData,N,wavename,DorS,parallelWorkers)
% Extract signal variations at a particular scale or set of scales using 
% MODWT (maximal overlap discrete wavelet transform).
%
% Usage: [waveData] = waveletTransform(rawData,N,wavename,DorS)
% 
% NOTE: The use of this code requires the WMTSA Toolbox available at: 
%       http://www.atmos.washington.edu/~wmtsa/
%
% ------------ Inputs -------------
% rawData = n x v data matrix to be decomposed, where n is the sample size
%      and v is the number of variables
% N = the dyadic (base-2) scale(s) to decompose.
%      For example, scale 5 identifies variations on the order of 2^5 data
%      points in width. N must be a single value unless detail
%      reconstruction is chosen (DorS = 1). In this case, the detail
%      reconstructions over the N scales will be summed.
% wavename = a string indicating the name of wavelet filter to use
%      Options are: 'haar', 'd4', 'la8', 'la16'
% DorS = choose to output (1) detail added at this scale (1), or (2) 
%      smooth approximation at this scale (low pass filter prior to 
%      adding detail at this scale)
% parallelWorkers = number of parallel workers to connect to if able to run 
%            script in parallel. Enter 1 to not run in parallel.
%
% ----------- Outputs ------------
% waveData = the wavelet transformed series, of same size as rawData;
%
% -------------------------------
% Cove Sturtevant, 2015.

plotfigs = 0; % plot figures while processing?

% Collect stats on data
sz = size(rawData);
if numel(sz) == 2
    sz(3) = 1;
end
nscales = length(N);

% Do a check on N
if nscales > 1 && DorS == 2
    N = N(1);
    logwrite(['Warning: wavelet scale must be single-valued for smooth approximation output. Using N = ' num2str(N)],1)
end
   
% Do check on NaNs
if sum(isnan(rawData(:))) ~= 0
    logwrite('Warning: Missing values found in data. Wavelet transform may be invalid.',1)
end
    
% Reset # of scales if it changed
nscales = length(N);

% Initialize output matrix     
waveData = NaN(sz(1),sz(2),sz(3));
    
% Perform wavelet transformation 
parfor (d2i = 1:sz(2),parallelWorkers)

    X = rawData(:,d2i);

    [DJt, SJt, mra_att] = modwt_mra(X, wavename, max(N),'reflection');

    % Remove reflected portion of series
    DJt = DJt(1:mra_att.NX,:);
    SJt = SJt(1:mra_att.NX,:);


    if plotfigs
        % Plot detail and approximation reconstructions
        figure(1); clf; 
        subplot(1,2,1); hold on
        yspread = (nanmax(X)-nanmin(X))*1.1;
        if yspread == 0
            yspread = 1;
        end
        yticks = nanmean(X):yspread:(nanmean(X)+yspread*nscales);
        Dplot = fliplr(DJt);
        for i = 1:nscales
            plot(Dplot(:,N(i))+yticks(i),'k')
        end
        axis tight
        set(gca,'ytick',yticks,'ylim',[yticks(1)-yspread/2 yticks(end)+yspread/2])
        set(gca,'yticklabel',[cellstr(num2str(flipud(N')));{'Data'}])
        title('Detail added at each scale')
        ylabel('Scale (lower scale = higher resolution)')

        subplot(1,2,2); hold on
        for i = 1:nscales+1
            if i == 1
                plot(SJt+yticks(i),'k')
            else
                plot(SJt+sum(Dplot(:,1:N(i-1)),2)+yticks(i),'k')
            end
        end
        axis tight
        set(gca,'ytick',yticks,'ylim',[yticks(1)-yspread/2 yticks(end)+yspread/2])
        set(gca,'yticklabel',[])
        title('Reconstructed time series at each scale')

        % Plot scale-wise variance contribution
        figure(2); clf; 
        bar(N,var(DJt(:,N),0,1))
        ylabel('Variance')
        xlabel('Scale')
        title(['Scale-wise variance partitioning: Var' num2str(d2i)])
        axis tight
        pause(0.5)

    end

    % Sum detail over scales, our output approximation
    if DorS == 1
        waveData(:,d2i) = sum(DJt(:,N),2);
    else
        waveData(:,d2i) = SJt;
    end
end
