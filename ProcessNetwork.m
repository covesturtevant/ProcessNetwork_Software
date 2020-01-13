function [R,opts] = ProcessNetwork(opts)
% 
% This code contains Version 1.5a of the main ProcessNetwork software, 
% coordinating all processing steps.
% 
% ------- Inputs -------
% opts = a structure with each field listing the parameters and settings
%        for processing. See Readme_ProcessNetwork_v1.5.docx for required
%        fields and defaults.
%
% ------- Outputs -------
% R = a structure of results. See Readme_ProcessNetwork_v1.5.docx for
%     descriptions.
% opts = the options structure input into the function, with defaults
%     filled into missing/bad fields.
% processLog = A cell of strings listing the progress of the processing 
%     steps performed. This is a global variable not output with the 
%     function, but can be accessed by typing the following from the
%     command window:
%     >> global processLog 
%
% -----------------------
% This version 1.5 of the main function was rewritted by Cove Sturtevant in
% July 2015, patterned after Version 1.4 written by Ben Ruddell. Vers. 1.5
% has a substantially different architecture compared to Version 1.4 to
% reduce computation time, add a wavelet tranformation module and IAAFT 
% surrogate data genereation, and repeat all processing steps similarly on
% the surrogate data.
% Version 1.5a adds a couple options and does further optimizing of the
% code for computation time
% - Added output for Linear Redundancy (Palus 2008) - See Sturtevant et al.
% 2016 for rationale
% - Added option to choose 'target' variables (akin to y variable in an x-y
% correlation). This is an option purely to reduce processing time by 
% not evaluating relationships for explanatory variables.
% - Added option to surrogate test each lag. This is important when there
% are gaps in the data.
% - Added output of relative mutual information normalized by (divided by) 
% the mean relative mutual information of the shuffled surrogates, and 
% associated statistical significance. New outputs are
% R.IRnormbyShuff
% R.SigThreshIRnormbyShuff
% R.sigmaShuffIRnormByShuff

%% Check parameters and settings
global processLog

clk = datestr(clock,'yyyy_mmm_dd_HH_MM_SS');
plogName = ['processLog_' clk]; % name of processLog for saving

logwrite(['Processing session for: ' clk],0);

logwrite('Checking options & parameters...',1);
opts = optsCheck(opts,[]);

% Intitialize
nDataFiles = size(opts.files,1);
badfile = 0; % skipped files
nLags = length(opts.lagVect);


%% Main program
logwrite('*** Beginning processing ***',1);

% Start the parallel pool
if opts.parallelWorkers > 1
    p = gcp;
    nWorkers = min(opts.parallelWorkers,p.NumWorkers);
    logwrite(['Matlab Pool open with ' num2str(nWorkers), ' CPU cores employed.'],1);
end

% Loop through files to process. 
for fi = 1:nDataFiles
    
    % Clear previously generated variables
    clear rawData Data Surrogates SavedSurrogates
        
    % Load file (can be matlab format or ascii
    logwrite(['--Processing file # ' num2str(fi) ': ' opts.files{fi} '...'],1);
    try
        [~,name,ext] = fileparts(opts.files{fi});
        if strcmp(ext,'.mat')
            load(opts.files{fi})
            size(Data);
        else
            Data = load(opts.files{fi});
        end
    catch
        logwrite('Error: Problem loading file. Skipping...',1);
        badfile = badfile+1;
        continue
    end
    [nData,nVars] = size(Data);
    
    % Retain data as loaded
    rawData = Data;
    
    % Check or create variable names, symbols, and units
    opts = optsCheck(opts,nVars);
    if length(opts.varNames) ~= nVars
        logwrite(['Unable to process file # ' num2str(fi) ': ' opts.files{fi} '. # of data columns inconsistent with # of varNames.'],1);
        badfile = badfile+1;
        continue
    end
    
    % Run the preprocessing options, including data trimming, anomaly
    % filter or wavelet transform
    logwrite('Preprocessing data.',1);
    if opts.trimTheData
        logwrite('Trimming rows with missing data',1);
    end
    if opts.transformation == 1
        logwrite(['Applying anomaly filter over ' num2str(opts.anomalyMovingAveragePeriodNumber) ' periods of ' num2str(opts.anomalyPeriodInData) ' time steps per period.'],1);
    elseif opts.transformation == 2
        if opts.waveDorS == 1
            DorS = 'detail';
        else
            DorS = 'approximation';
        end
        logwrite(['Applying MODWT wavelet filter at ' DorS ' scale(s) [' num2str(opts.waveN) '] using ' opts.waveName ' mother wavelet.'],1);
    end
    Data = preProcess(opts,Data);

    % Intialize if first file
    if fi-badfile == 1
        R = intializeOutput(nDataFiles,nVars,opts);
    end
    
    % Collect summary stats on variables for archiving and later
    % classification
    logwrite('Computing local statistics',1);
    R.nRawData(fi) = nData;
    R.nVars(fi) = nVars;
    [R.binEdgesLocal(:,:,fi),R.minEdgeLocal(:,fi),R.maxEdgeLocal(:,fi)] = GetUniformBinEdges(Data,R.nBinVect,opts.binPctlRange,NaN);
    R.LocalVarAvg(:,fi) = nanmean(Data,1)';
    R.LocalVarCnt(:,fi) = sum(~isnan(Data))';

    % If we are doing local binning or data is already binned, we can 
    % go straight into classification and/or entropy calculations. 
    % Otherwise, we are just saving stats for now
    if ~isempty(find([0 1] == opts.binType,1))
        
        % Classify the data with local binning
        if opts.binType == 1
            logwrite(['Classifying with [' num2str(opts.nBins) '] local bins over [' num2str(opts.binPctlRange) '] percentile range.'],1);
            [Data,R.nClassified(fi)]=classifySignal(Data,R.binEdgesLocal(:,:,fi),R.nBinVect,NaN);
        end
        
        % Save preprocessed Data if we're done processing it
        if opts.savePreProcessed && ~isempty(find([0 1] == opts.binType,1))
            logwrite('Saving preprocessed data',1);
            eval([plogName ' = processLog;'])
            save([opts.outDirectory name opts.preProcessedSuffix],'Data',plogName)
        end
        
        % Run entropy function
        if opts.doEntropy
            logwrite('Running entropy function.',1);
            
            % Check that data have been classified
            if nansum(rem(Data(:),1)) ~= 0
                logwrite('ERROR: Some or all input data are not classified. Check options and/or data. Skipping file...',1);
                continue
            end
            
            [E] = entropyFunction(Data,opts.targetVars,R.lagVect,R.nBinVect,NaN,opts.parallelWorkers);

            % Assign outputs
            R.HXt(:,:,:,fi)=E.HXt;
            R.HYw(:,:,:,fi)=E.HYw;
            R.HYf(:,:,:,fi)=E.HYf;
            R.HXtYw(:,:,:,fi)=E.HXtYw;
            R.HXtYf(:,:,:,fi)=E.HXtYf;
            R.HYwYf(:,:,:,fi)=E.HYwYf;
            R.HXtYwYf(:,:,:,fi)=E.HXtYwYf;
            R.nCounts(:,:,:,fi)=E.nCounts;
            R.I(:,:,:,fi)=E.I;
            R.T(:,:,:,fi)=E.T;
            R.L(:,:,:,fi)=E.L;
        end
    end
    
    % If we already have surrogates, check we have them
    if opts.SurrogateMethod == 1
        if ~exist('Surrogates','var')
            msg = logwrite('FATAL ERROR: No Surrogates exist in loaded file. Check SurrogateMethod option.',0);
            error(msg)
        elseif size(Surrogates,3) < opts.nTests
            msg = logwrite('FATAL ERROR: Number of loaded Surrogates is less than nTests in options.',0);
            error(msg)
        end
        SavedSurrogates = Surrogates;
    elseif opts.savePreProcessed == 1 && opts.SurrogateMethod > 1
        % Initialize if we are saving them and don't already have them
        SavedSurrogates = NaN([nData nVars opts.nTests]);
    end
        
    % Create and/or process surrogates for statistical significance 
    % testing
    if opts.SurrogateMethod > 0
        
        if opts.SurrogateMethod == 1
            logwrite(['Running the same operations on ' num2str(opts.nTests) ' surrogates contained in input file...'],1);
        elseif opts.SurrogateMethod == 2
            logwrite(['Creating and running the same operations on ' num2str(opts.nTests) ' surrogates using random shuffle method.'],1);
        elseif opts.SurrogateMethod == 3
            logwrite(['Creating and running the same operations on ' num2str(opts.nTests) ' surrogates using IAAFT method (this may take a while).'],1);
        elseif opts.SurrogateMethod == 4
            logwrite(['Creating and running the same operations on ' num2str(opts.nTests) ' surrogates using random walk method (this may take a while).'],1);
        end
        
        % Initalize surrogate matrix
        if opts.SurrogateTestEachLag == 1
            % Test every lag in surrogates
            SlagVect = R.lagVect;
            nSLags = nLags;
            if ~isempty(find([0 1] == opts.binType,1)) && opts.doEntropy == 1
                logwrite('Testing surrogates at each lag (this may take a while).',1);
            end
        else
            % Test only the first and last lags (to restrict data range
            % to same as actual data). We will only retain stats for the
            % last lag.
            SlagVect = [R.lagVect(1) R.lagVect(end)]';
            nSLags = 1;
            if ~isempty(find([0 1] == opts.binType,1)) && opts.doEntropy == 1
                logwrite('Testing surrogates at final lag only.',1);
            end
        end
        shuffT = NaN([nVars nVars nSLags opts.nTests]);
        shuffI = NaN([nVars nVars nSLags opts.nTests]);
        shuffL = NaN([nVars nVars nSLags opts.nTests]);
        shuffHYf = NaN([nVars nVars nSLags opts.nTests]);        

        for si = 1:opts.nTests
            if opts.SurrogateMethod == 1
                Surrogates = SavedSurrogates(:,:,si);
            elseif ~isempty(find([2 3 4] == opts.SurrogateMethod,1))
                % Create surrogates using method specified
                Surrogates = createSurrogates(opts,rawData,1);
            end
            
            % Preprocess surrogates same as Data
            Surrogates = preProcess(opts,Surrogates);
            
            % Collect stats on Surrogates
            [SbinEdgesLocal,minEdgeLocal,maxEdgeLocal] = GetUniformBinEdges(Surrogates,R.nBinVect,opts.binPctlRange,NaN);
            R.minSurrEdgeLocal(:,fi) = nanmin([minEdgeLocal R.minSurrEdgeLocal(:,fi)],[],2);
            R.maxSurrEdgeLocal(:,fi) = nanmax([maxEdgeLocal R.maxSurrEdgeLocal(:,fi)],[],2);
            
            % If we are doing local binning or data is already binned, we can 
            % go straight into classification and/or entropy calculations. 
            % Otherwise, we are just saving stats for now
            if ~isempty(find([0 1] == opts.binType,1))

                % Classify the data with local binning
                if opts.binType == 1
                    [Surrogates,~]=classifySignal(Surrogates,SbinEdgesLocal,R.nBinVect,NaN);
                end

                % Are we saving the surrogates?
                if opts.savePreProcessed
                    SavedSurrogates(:,:,si) = Surrogates;
                end
                
                % Run entropy function
                if opts.doEntropy
                    
                    if nansum(rem(Surrogates(:),1)) ~= 0
                        logwrite('ERROR: Surrogate data are not classified. Check options and/or input data. Aborting surrogate testing...',1);
                        break
                    end
                    [E] = entropyFunction(Surrogates,opts.targetVars,SlagVect,R.nBinVect,NaN,opts.parallelWorkers);
                    
                    % Assign outputs specific to surrogate data
                    if opts.SurrogateTestEachLag == 1
                        % All lags tested
                        shuffT(:,:,:,si) = E.T;
                        shuffI(:,:,:,si) = E.I;
                        shuffL(:,:,:,si) = E.L;
                        shuffHYf(:,:,:,si) = E.HYf;
                        
                    else
                        % Just last lag
                        shuffT(:,:,1,si) = E.T(:,:,end);
                        shuffI(:,:,1,si) = E.I(:,:,end);
                        shuffL(:,:,1,si) = E.L(:,:,end);
                        shuffHYf(:,:,1,si) = E.HYf(:,:,end);
                    end
                end
            end
        end
        
        % Calculate stats for statistical significance
        R.meanShuffT(:,:,:,fi)=mean(shuffT,4);
        R.sigmaShuffT(:,:,:,fi)=std(shuffT,0,4);
        R.meanShuffI(:,:,:,fi)=mean(shuffI,4);
        R.sigmaShuffI(:,:,:,fi)=std(shuffI,0,4);
        R.meanShuffL(:,:,:,fi)=mean(shuffL,4);
        R.sigmaShuffL(:,:,:,fi)=std(shuffL,0,4);
        R.meanShuffTR(:,:,:,fi)=mean(shuffT./shuffHYf,4);
        R.sigmaShuffTR(:,:,:,fi)=std(shuffT./shuffHYf,0,4);
        R.meanShuffIR(:,:,:,fi)=mean(shuffI./shuffHYf,4);
        R.sigmaShuffIR(:,:,:,fi)=std(shuffI./shuffHYf,0,4);

       % Save preprocessed Surrogates if we're done processing
        if opts.savePreProcessed && ~isempty(find([0 1] == opts.binType,1))
            logwrite('Saving preprocessed surrogates',1);
            Surrogates = SavedSurrogates;
            eval([plogName ' = processLog;'])
            save([opts.outDirectory name opts.preProcessedSuffix],'Surrogates',plogName,'-append')
        end
    end
end
                
%Make sure we processed at least 1 file
if badfile == nDataFiles
    logwrite('No files were processed. Check processLog',1);
    R = [];
    return
end

% Establish global statistics and bins
logwrite('Computing global statistics',1);
R.GlobalVarAvg = nansum(R.LocalVarAvg.*R.LocalVarCnt,2)./nansum(R.LocalVarCnt,2);
[R.binEdgesGlobal,R.minEdgeGlobal,R.maxEdgeGlobal]=GetEvenBinEdgesGlobal(R.nBinVect,R.minEdgeLocal,R.maxEdgeLocal); % data
[R.binSurrEdgesGlobal,R.minSurrEdgeGlobal,R.maxSurrEdgeGlobal]=GetEvenBinEdgesGlobal(R.nBinVect,nanmin(R.minSurrEdgeLocal,[],3),nanmax(R.maxSurrEdgeLocal,[],3)); % surrogates


% If we chose the global binning option, we need to run through the data
% again
if opts.binType == 2
    logwrite('*** Processing files again, this time using global binning ***',1);

    % Reset Surrogate stats
    R.minSurrEdgeLocal = NaN(nVars,nDataFiles);
    R.maxSurrEdgeLocal = NaN(nVars,nDataFiles);

    for fi = 1:nDataFiles
        logwrite(['--Processing file # ' num2str(fi) ': ' opts.files{fi} '...'],1);

        % Clear previously generated variables
        clear rawData Data Surrogates SavedSurrogates

        % Load file (can be matlab format or ascii
        try
            [~,name,ext] = fileparts(opts.files{fi});
            if strcmp(ext,'.mat')
                load(opts.files{fi})
                if ~exist('Data','var')
                    error('')
                end
            else
                Data = load(opts.files{fi});
            end
        catch
            logwrite(['Unable to process file # ' num2str(fi) ': ' opts.files{fi} '. Problem loading file.'],1)
            badfile = badfile+1;
            continue
        end
        [nData,nVars] = size(Data);

        % Retain data as loaded
        rawData = Data;

        % Run the preprocessing options, including data trimming, anomaly
        % filter or wavelet transform
        logwrite('Preprocessing data.',1);
        if opts.trimTheData
            logwrite('Trimming rows with missing data',1);
        end
        if opts.transformation == 1
            logwrite(['Applying anomaly filter over ' num2str(opts.anomalyMovingAveragePeriodNumber) ' periods of ' num2str(opts.anomalyPeriodInData) ' time steps per period.'],1);
        elseif opts.transformation == 2
            if opts.waveDorS == 1
                DorS = 'detail';
            else
                DorS = 'approximation';
            end
            logwrite(['Applying MODWT wavelet filter at ' DorS ' scale(s) [' num2str(opts.waveN) '] using ' opts.waveName ' mother wavelet.'],1);
        end
        Data = preProcess(opts,Data);

        % Classify the data with global binning
        logwrite(['Classifying with [' num2str(opts.nBins) '] global bins over [' num2str(opts.binPctlRange) '] percentile range.'],1);
        [Data,R.nClassified(fi)]=classifySignal(Data,R.binEdgesGlobal,R.nBinVect,NaN);

        % Save preprocessed Data
        if opts.savePreProcessed
            logwrite('Saving preprocessed data',1);
            eval([plogName ' = processLog;'])
            save([opts.outDirectory name opts.preProcessedSuffix],'Data',plogName)
        end
        
        % Run entropy function
        if opts.doEntropy
            logwrite('Running entropy function.',1);
            [E] = entropyFunction(Data,opts.targetVars,R.lagVect,R.nBinVect,NaN,opts.parallelWorkers);

            % Assign outputs
            R.HXt(:,:,:,fi)=E.HXt;
            R.HYw(:,:,:,fi)=E.HYw;
            R.HYf(:,:,:,fi)=E.HYf;
            R.HXtYw(:,:,:,fi)=E.HXtYw;
            R.HXtYf(:,:,:,fi)=E.HXtYf;
            R.HYwYf(:,:,:,fi)=E.HYwYf;
            R.HXtYwYf(:,:,:,fi)=E.HXtYwYf;
            R.nCounts(:,:,:,fi)=E.nCounts;
            R.I(:,:,:,fi)=E.I;
            R.T(:,:,:,fi)=E.T;
            R.L(:,:,:,fi)=E.L;
        end

        % Are we saving the Surrogates? If so, load or initialize matrix
        if opts.SurrogateMethod == 1
            SavedSurrogates = Surrogates;
        elseif opts.savePreProcessed == 1 && opts.SurrogateMethod > 1
            SavedSurrogates = NaN([nData nVars opts.nTests]);
        end

        % Create and/or process surrogates for statistical significance 
        % testing
        if opts.SurrogateMethod > 0
        
            if opts.SurrogateMethod == 1
                logwrite(['Running the same operations on ' num2str(opts.nTests) ' surrogates contained in input file (this may take a while)...'],1);
            elseif opts.SurrogateMethod == 2
                logwrite(['Creating and running the same operations on ' num2str(opts.nTests) ' surrogates using random shuffle method (this may take a while)...'],1);
            elseif opts.SurrogateMethod == 3
                logwrite(['Creating and running the same operations on ' num2str(opts.nTests) ' surrogates using IAAFT method (this may take a while)...'],1);
            elseif opts.SurrogateMethod == 4
                logwrite(['Creating and running the same operations on ' num2str(opts.nTests) ' surrogates using random walk method (this may take a while)...'],1);
            end

            % Initalize surrogate matrix
            if opts.SurrogateTestEachLag == 1 && opts.doEntropy == 1
                % Test every lag in surrogates
                SlagVect = R.lagVect;
                nSLags = nLags;
                logwrite('Testing surrogates at each lag (this may take a while).',1);
            elseif opts.doEntropy == 1
                % Test only the first and last lags (to restrict data range
                % to same as actual data). We will only retain stats for the
                % last lag.
                SlagVect = [R.lagVect(1) R.lagVect(end)]';
                nSLags = 1;
                logwrite('Testing surrogates at final lag only.',1);
            end
            shuffT = NaN([nVars nVars nSLags opts.nTests]);
            shuffI = NaN([nVars nVars nSLags opts.nTests]);
            shuffL = NaN([nVars nVars nSLags opts.nTests]);
            shuffHYf = NaN([nVars nVars nSLags opts.nTests]);        
        
            for si = 1:opts.nTests
                if opts.SurrogateMethod == 1
                    Surrogates = SavedSurrogates(:,:,si);
                elseif ~isempty(find([2 3 4] == opts.SurrogateMethod,1))
                    % Create surrogates using method specified
                    [Surrogates] = createSurrogates(opts,rawData,1);
                end

                % Preprocess surrogates same as Data
                Surrogates = preProcess(opts,Surrogates);

                % Collect stats on Surrogates - we're going to rewrite
                % these to ensure that the global stats we calculated from
                % different surrogates match the new ones we create here
                [~,minEdgeLocal,maxEdgeLocal] = GetUniformBinEdges(Surrogates,R.nBinVect,opts.binPctlRange,NaN);
                R.minSurrEdgeLocal(:,fi) = nanmin([minEdgeLocal R.minSurrEdgeLocal(:,fi)],[],2);
                R.maxSurrEdgeLocal(:,fi) = nanmax([maxEdgeLocal R.maxSurrEdgeLocal(:,fi)],[],2);

                % Classify surrogates
                [Surrogates,~]=classifySignal(Surrogates,R.binSurrEdgesGlobal,R.nBinVect,NaN);

                % Are we saving the surrogates? If so, archive.
                if opts.savePreProcessed
                    SavedSurrogates(:,:,si) = Surrogates;
                end

                % Run entropy function
                if opts.doEntropy
                    
                    [E] = entropyFunction(Surrogates,opts.targetVars,SlagVect,R.nBinVect,NaN,opts.parallelWorkers);
                    
                    % Assign outputs specific to surrogate data
                    if opts.SurrogateTestEachLag == 1
                        % All lags tested
                        shuffT(:,:,:,si) = E.T;
                        shuffI(:,:,:,si) = E.I;
                        shuffL(:,:,:,si) = E.L;
                        shuffHYf(:,:,:,si) = E.HYf;
                        
                    else
                        % Just last lag
                        shuffT(:,:,1,si) = E.T(:,:,end);
                        shuffI(:,:,1,si) = E.I(:,:,end);
                        shuffL(:,:,1,si) = E.L(:,:,end);
                        shuffHYf(:,:,1,si) = E.HYf(:,:,end);
                    end
                end              
            end

            % Calculate stats for statistical significance
            R.meanShuffT(:,:,:,fi)=mean(shuffT,4);
            R.sigmaShuffT(:,:,:,fi)=std(shuffT,0,4);
            R.meanShuffI(:,:,:,fi)=mean(shuffI,4);
            R.sigmaShuffI(:,:,:,fi)=std(shuffI,0,4);
            R.meanShuffL(:,:,:,fi)=mean(shuffL,4);
            R.sigmaShuffL(:,:,:,fi)=std(shuffL,0,4);
            R.meanShuffTR(:,:,:,fi)=mean(shuffT./shuffHYf,4);
            R.sigmaShuffTR(:,:,:,fi)=std(shuffT./shuffHYf,0,4);
            R.meanShuffIR(:,:,:,fi)=mean(shuffI./shuffHYf,4);
            R.sigmaShuffIR(:,:,:,fi)=std(shuffI./shuffHYf,0,4);

            % Save preprocessed Surrogates
            if opts.savePreProcessed
                logwrite('Saving preprocessed surrogates',1);
                Surrogates = SavedSurrogates;
                eval([plogName ' = processLog;'])
                save([opts.outDirectory name opts.preProcessedSuffix],'Surrogates',plogName,'-append')
            end
        end
    end
end

% Close the parallel pool
if opts.parallelWorkers > 1
    if opts.closeParallelPool
        delete(gcp('nocreate'))
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% POSTPROCESS DERIVED INFORMATION THEORY AND PHYSICAL QUANTITIES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Statistical significance thresholds
if opts.SurrogateMethod > 0
    % Statistical significance thresholds
    R.SigThreshT = R.meanShuffT+opts.oneTailZ*R.sigmaShuffT;
    R.SigThreshI = R.meanShuffI+opts.oneTailZ*R.sigmaShuffI;
    R.SigThreshL = R.meanShuffL+opts.oneTailZ*R.sigmaShuffL;
    R.SigThreshTR = R.meanShuffTR+opts.oneTailZ*R.sigmaShuffTR;
    R.SigThreshIR = R.meanShuffIR+opts.oneTailZ*R.sigmaShuffIR;
    
end

% Derived Quantities
if opts.doEntropy
    logwrite('Computing final entropy quantities.',1);
    [R.Tplus,R.Tminus,R.Tnet,R.TnetBinary]=DoProduction(R.T);
    [R.InormByDist,R.TnormByDist,R.SigThreshInormByDist,R.SigThreshTnormByDist,R.Ic,R.Tc,R.TvsIzero,R.SigThreshTvsIzero,R.IR,R.TR,R.HXtNormByDist,R.IvsIzero,R.SigThreshIvsIzero]=NormTheStats(R.nBinVect,R.I,R.T,R.SigThreshI,R.SigThreshT,R.meanShuffI,R.sigmaShuffI,R.meanShuffT,R.sigmaShuffT,R.HXt,R.HYf,R.lagVect);
    [R.Abinary,R.Awtd,R.AwtdCut,R.charLagFirstPeak,R.TcharLagFirstPeak,R.charLagMaxPeak,R.TcharLagMaxPeak,R.TvsIzerocharLagMaxPeak,R.nSigLags,R.FirstSigLag,R.LastSigLag]=AdjMatrices(R.TnormByDist,R.SigThreshTnormByDist,R.TvsIzero,R.lagVect);
    R.Hm=sum(squeeze(R.HXt(:,1,1,:)./repmat(log2(R.nBinVect),[1 1,1,nDataFiles])))./nVars;
    R.TSTm=squeeze(sum(sum(R.T./repmat(log2(R.nBinVect),[1 nVars nLags nDataFiles]),1),2))./(nVars^2);
    
    % Normalization with mean IR of surrogates
    if opts.SurrogateTestEachLag == 1
        R.IRnormByShuff = R.IR./R.meanShuffIR;
        R.sigmaShuffIRnormByShuff = R.sigmaShuffIR./R.meanShuffIR;
        R.SigThreshIRnormByShuff = 1+opts.oneTailZ*R.sigmaShuffIRnormByShuff;
    end

end       

% Save output
if opts.saveProcessNetwork == 1
    logwrite('Saving results.',1);
    if isnan(opts.outFileProcessNetwork)
        outfile = ['Results_ ' clk];
    else
        outfile = opts.outFileProcessNetwork;
    end
    
    % Rename the processLog so it is unique
    eval([plogName ' = processLog;'])
    save([opts.outDirectory outfile],'R','opts',plogName)
end

logwrite('Processing run complete.',1);

