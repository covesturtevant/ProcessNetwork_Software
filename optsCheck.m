function [opts] = optsCheck(opts,nVars)
% Checks parameters and settings to run the ProcessNetwork software
%
% ------- Inputs -------
% opts = the options structure listing parameters and settings (see README)
% nVars = The number of variables (columns) in the input data. If unknown, 
%         leave as NaN or empty []
% 
% ------- Output -------
% opts = the validated opts structure updated with default values where 
%        possible
% 
% ----------------------
if isempty(nVars)
    nVars = NaN;
end

if ~isfield(opts,'files')
    msg = logwrite('FATAL ERROR: No files listed to process. List files to process in options structure.',0);
    error(msg)
end

if ~isfield(opts,'varNames')
    if ~isnan(nVars)
        logwrite('Setting variable names to (V1,V2,etc.) for rest of processing run.',1);
        opts.varNames = cell([1 nVars]);
        for i = 1:nVars
            opts.varNames{i} = ['V' num2str(i)];
        end
    else
        logwrite('No varNames option found. Default names will be used (V1,V2,etc.)',1);
    end
elseif isempty(opts.varNames)
    if ~isnan(nVars)
        logwrite('Setting variable names to (V1,V2,etc.) for rest of processing run.',1);
        opts.varNames = cell([1 nVars]);
        for i = 1:nVars
            opts.varNames{i} = ['V' num2str(i)];
        end
    else
        logwrite('No varNames option found. Default names will be used (V1,V2,etc.)',1);
    end    
elseif ~iscell(opts.varNames)
    if ~isnan(nVars)
        logwrite('Setting variable names to (V1,V2,etc.) for rest of processing run.',1);
        opts.varNames = cell([1 nVars]);
        for i = 1:nVars
            opts.varNames{i} = ['V' num2str(i)];
        end
    else
        logwrite('Bad varNames option. Default names will be used (V1,V2,etc.)',1);
        opts = rmfield(opts,'varNames');
    end
end

if ~isfield(opts,'targetVars')
    if ~isnan(nVars)
        logwrite('Setting targetVars to all.',1);
        opts.targetVars = NaN([1 nVars]);
        for i = 1:nVars
            opts.targetVars(i) = i;
        end
    else
        logwrite('No targetVars option found. Setting targetVars to all.',1);
    end
elseif isempty(opts.targetVars)
    if ~isnan(nVars)
        logwrite('Setting targetVars to all for rest of processing run.',1);
        opts.targetVars = NaN([1 nVars]);
        for i = 1:nVars
            opts.targetVars(i) = i;
        end
    else
        logwrite('No targetVars option found. Setting targetVars to all for rest of processing run.',1);
    end    
elseif iscell(opts.targetVars) && ischar(cell2mat(opts.targetVars))
    targetVars = cell([1 length(opts.targetVars(:))]);
    for i=1:length(targetVars)
        targetVars{i} = find(strcmp(opts.varNames,opts.targetVars{i}));
    end
    opts.targetVars=targetVars;
    
    if length(cell2mat(opts.targetVars)) ~= length(opts.targetVars)
        opts.targetVars = NaN([1 nVars]);
        for i = 1:nVars
            opts.targetVars(i) = i;
        end
        logwrite('At least one variable listed in targetVars does not match the list of varNames. Setting targetVars to all.',1);
    else
        opts.targetVars=unique(cell2mat(opts.targetVars));
    end
elseif iscell(opts.targetVars) && ~isnan(nVars) && max(cell2mat(opts.targetVars(:))) > nVars
     opts.targetVars = NaN([1 nVars]);
        for i = 1:nVars
            opts.targetVars(i) = i;
        end
    logwrite('At least one indices in targetVars is greater than the number of variables. Setting targetVars to all.',1);
elseif ~iscell(opts.targetVars) && ~isnan(nVars) && max(opts.targetVars(:)) > nVars
     opts.targetVars = NaN([1 nVars]);
        for i = 1:nVars
            opts.targetVars(i) = i;
        end
    logwrite('At least one indices in targetVars is greater than the number of variables. Setting targetVars to all.',1);
end

if ~isfield(opts,'varSymbols')
    if ~isnan(nVars)
        logwrite('Setting variable symbols to variable names for rest of processing run.',1);
        opts.varSymbols = opts.varNames;
    else
        logwrite('No varSymbols option found. varNames will be used in labeling.',1);
    end
elseif (length(opts.varSymbols) ~= nVars || ~iscell(opts.varSymbols) )&& (~isnan(nVars))
    logwrite('Warning: # of varSymbols inconsistent with # of varNames. Setting varSymbols to varNames.',1);
    opts.varSymbols = opts.varNames;    
end

if ~isfield(opts,'varUnits')
    if ~isnan(nVars)
        opts.varUnits = cell([1 nVars]);
    else
        logwrite('No varUnits option found. No variable units will be displayed.',1);
    end
elseif (length(opts.varUnits) ~= nVars || ~iscell(opts.varUnits)) && (~isnan(nVars))
    logwrite('Warning: # of varUnits inconsistent with # of varNames. Clearing varUnits.',1);
    opts.varUnits = cell([1 nVars]);
end

if ~isfield(opts,'NoDataCode')
    logwrite('No NoDataCode option found. Setting to default (NaN).',1);
    opts.NoDataCode = NaN; % data with this value will be skipped in entropy computations. Must be NaN or an integer that can't possibly be a bin #, i.e. negative integer! default  = NaN.
elseif ~isnumeric(opts.NoDataCode)
    logwrite('Warning: Bad NoDataCode option. Setting to default (NaN).',1);
    opts.NoDataCode = NaN; % data with this value will be skipped in entropy computations. Must be NaN or an integer that can't possibly be a bin #, i.e. negative integer! default  = NaN.
end

if ~isfield(opts,'trimTheData') 
    logwrite('No trimTheData option found. Setting to default (0 = no).',1);
    opts.trimTheData = 0; % set to 0 to not trim the whole row if one data value in the row is NoData; default is 1
elseif isempty(find([0 1] == opts.trimTheData,1))
    logwrite('Warning: Bad trimTheData option. Setting to default (0 = no).',1);
    opts.trimTheData = 0; % set to 0 to not trim the whole row if one data value in the row is NoData; default is 1
end

if ~isfield(opts,'transformation')
    logwrite('No transformation option found. Setting to default (0 = none).',1);
    opts.transformation = 0; %set to 0 to not apply any tranformation, 1 for anomaly filter, 2 for wavelet transform 
elseif isempty(find([0 1 2] == opts.transformation,1))
    logwrite('Warning: Bad transformation option. Setting to default (0 = none).',1);
    opts.transformation = 0; %set to 0 to not apply any tranformation, 1 for anomaly filter, 2 for wavelet transform 
end

% Anomaly filter settings
if opts.transformation == 1
    if ~isfield(opts,'anomalyPeriodInData')
        logwrite('No anomalyPeriodInData option found. Setting to default (48).',1);
        opts.anomalyPeriodInData = 48; %set to the length in time steps of the period of the data
    elseif ~isnumeric(opts.anomalyPeriodInData)
        logwrite('Warning: Bad anomalyPeriodInData option. Setting to default (48).',1);
        opts.anomalyPeriodInData = 48; %set to the length in time steps of the period of the data
    end

    if ~isfield(opts,'anomalyMovingAveragePeriodNumber')
        logwrite('No anomalyMovingAveragePeriodNumber option. Setting to default (5).',1);
        opts.anomalyMovingAveragePeriodNumber = 5; % the moving window used for anomaly generation, in #periods
    elseif ~isnumeric(opts.anomalyMovingAveragePeriodNumber)
        logwrite('Warning: Bad anomalyMovingAveragePeriodNumber option. Setting to default (5).',1);
        opts.anomalyMovingAveragePeriodNumber = 5; % the moving window used for anomaly generation, in #periods
    end
end

% Wavelet filter settings
if opts.transformation == 2
    if ~isfield(opts,'waveN') && opts.transformation == 2
        logwrite('No waveN option found. Setting to default (1).',1);
        opts.waveN = 1; % 1 x n numeric array of dyadic wavelet scales to decompose
    elseif opts.transformation == 2 && ~isnumeric(opts.waveN)
        logwrite('Warning: Bad waveN option. Setting to default (1).',1);
        opts.waveN = 1; % 1 x n numeric array of dyadic wavelet scales to decompose
    end

    if ~isfield(opts,'waveName')
        logwrite('No waveName option found. Setting to default (la8).',1);
        opts.waveName = 'la8'; % string indicating mother wavelet to use in wavelet transform (defailt = 'la8')
    end
    if opts.transformation == 2 && ~strcmp(opts.waveName,'haar') && ~strcmp(opts.waveName,'d4') && ~strcmp(opts.waveName,'la8') && ~strcmp(opts.waveName,'la16')
        logwrite('Warning: waveName not recognized. Setting to default (la8).',1);
        opts.waveName = 'la8'; % string indicating mother wavelet to use in wavelet transform (defailt = 'la8')
    end

    if ~isfield(opts,'waveDorS')
        logwrite('No waveDorS option found. Setting to default (1 = detail).',1);
        opts.waveDorS = 1; %output (1) detail added at this scale (1), or (2) approximation
    elseif isempty(find([1 2] == opts.waveDorS,1))
        logwrite('Warning: Bad waveDorS option. Setting to default (1 = detail).',1);
        opts.waveDorS = 1; %output (1) detail added at this scale (1), or (2) approximation
    end
    if opts.transformation == 2 && length(opts.waveN) > 1
        if opts.waveDorS == 2
            opts.waveN = opts.waveN(1);
            logwrite(['Warning: waveN must be single-valued for smooth approximation output. Truncating waveN to ' num2str(opts.waveN)],1);
        end
    end
end

if ~isfield(opts,'binType')
    logwrite('No binType option found. Setting to default (1 = do local binning).',1);
    opts.binType = 1; % data binning: 0 = don't do classification (or data are already binned), 1 = classify using locally bounded bins, 2  = classify using globally bounded bins
elseif isempty(find([0 1 2] == opts.binType,1))
    logwrite('Warning: Bad binType option. Setting to default (1 = do local binning).',1);
    opts.binType = 1; % data binning: 0 = don't do classification (or data are already binned), 1 = classify using locally bounded bins, 2  = classify using globally bounded bins
end

% Classification
if opts.binType > 0
    if ~isfield(opts,'binPctlRange')
        logwrite('No binPctlRange option found. Setting to default ([0 100]).',1);
        opts.binPctlRange = [0 100]; % [min max] percentile range of data used to determine total bin range for each variable. Default is entire data range ([0 100]) 
    elseif numel(opts.binPctlRange) ~= 2 || min(opts.binPctlRange) < 0 || max(opts.binPctlRange) > 100
        logwrite('Warning: Bad binPctlRange option. Setting to default ([0 100]).',1);
        opts.binPctlRange = [0 100]; % [min max] percentile range of data used to determine total bin range for each variable. Default is entire data range ([0 100]) 
    end
end

if ~isfield(opts,'nBins')
    logwrite('No nBins option found. Setting to default (11).',1);
    opts.nBins = 11; % set the number of bins used to classify each signal here. 
elseif ~isnan(nVars) && length(opts.nBins) < nVars && length(opts.nBins) > 1
    logwrite('Length of nBins does not match number of variables. Setting to default (11).',1);
    opts.nBins = 11; % set the number of bins used to classify each signal here. 
elseif ~isnumeric(opts.nBins) || min(opts.nBins) < 2
    logwrite('Warning: Bad nBins option. Setting to default (11).',1);
    opts.nBins(opts.nBins < 2) = 11; % set the number of bins used to classify each signal here
end
if size(opts.nBins,1) > size(opts.nBins,2)
    opts.nBins = opts.nBins'; % Make row vector
end

if ~isfield(opts,'SurrogateMethod')
    logwrite('No SurrogateMethod option found. Setting to default (2 = random shuffle).',1);
    opts.SurrogateMethod = 2; % 0 = No statistical testing, regardless of whether input files contain Surrogates. 1 = Use the Surrogates contained in the loaded files. 2 = new surrogates created via random shuffle of input data; 3 = new surrogates created via IAAFT method performed on input data
elseif isempty(find([0 1 2 3 4 5] == opts.SurrogateMethod,1))
    logwrite('Warning: Bad SurrogateMethod option. Setting to default (2 = random shuffle).',1);
    opts.SurrogateMethod = 2; % 0 = No statistical testing, regardless of whether input files contain Surrogates. 1 = Use the Surrogates contained in the loaded files. 2 = new surrogates created via random shuffle of input data; 3 = new surrogates created via IAAFT method performed on input data
end

if opts.SurrogateMethod > 0
    if ~isfield(opts,'nTests')
        logwrite('No nTests option found. Setting to default (100).',1);
        opts.nTests = 100; % Number of surrogate data series to test for determining statistical significance thresholds
    elseif ~isnumeric(opts.nTests)
        logwrite('Warning: Bad nTests option. Setting to default (100).',1);
        opts.nTests = 100; % Number of surrogate data series to test for determining statistical significance thresholds
    end
    
    if ~isfield(opts,'SurrogateTestEachLag')
        logwrite('No SurrogateTestEachLag option found. Setting to default (0 = no, tests last lag only).',1);
        opts.SurrogateTestEachLag = 0; % Provide statsistical significance threshold for each lag? 0 = no, test first lag only, 1 = yes, test each lag 
    elseif isempty(find([0 1] == opts.SurrogateTestEachLag,1))
        logwrite('Warning: Bad SurrogateTestEachLag option. Setting to default (0 = no, tests last lag only).',1);
        opts.SurrogateTestEachLag = 0; % Provide statsistical significance threshold for each lag? 0 = no, test first lag only, 1 = yes, test each lag 
    end
    
    if ~isfield(opts,'oneTailZ')
        logwrite('No oneTailZ option found. Setting to default (1.66).',1);
        opts.oneTailZ = 1.66; % one-tail z-score for 95% given number of tests, 1.66 for 100, 1.68 for 50, 1.71 for 25
    elseif ~isnumeric(opts.oneTailZ) || opts.oneTailZ < 0
        logwrite('Warning: Bad oneTailZ option. Setting to default (1.66).',1);
        opts.oneTailZ = 1.66; % one-tail z-score for 95% given number of tests, 1.66 for 100, 1.68 for 50, 1.71 for 25
    end
end

if ~isfield(opts,'doEntropy')
    logwrite('No doEntropy option found. Setting to default (0 = no).',1);
    opts.doEntropy = 0; % run entropy calculations?
elseif isempty(find([0 1] == opts.doEntropy,1))
    logwrite('Warning: Bad doEntropy option. Setting to default (0 = no).',1);
    opts.doEntropy = 0; % run entropy calculations?
end

% Lag checking
if opts.doEntropy == 1
    if ~isfield(opts,'lagVect')
        logwrite('No lagVect option found. Setting to default (0:10).',1);
        opts.lagVect = (0:10)'; % lags to evaluate in units of data time step
    elseif ~isnumeric(opts.lagVect)
        logwrite('Warning: Bad lagVect option. Setting to default (0:10).',1);
        opts.lagVect = (0:10)'; % lags to evaluate in units of data time step
    elseif isempty(find(opts.lagVect == 0,1))
        opts.lagVect(end+1) = 0;
    end
    
    if size(opts.lagVect,2) > size(opts.lagVect,1)
        opts.lagVect = opts.lagVect';
    end
    
    if size(opts.lagVect,2) > 1
        opts.lagVect = opts.lagVect(:,1);
        logwrite('Warning: lagVect option should be a single column vector. Truncating to first column.',1);
    end
    
    % Sort lags
    opts.lagVect = sort(opts.lagVect);
else
    opts.lagVect = 0;
end

if ~isfield(opts,'parallelWorkers')
    logwrite('No parallelWorkers option found. Setting to default (1 = no parallelization).',1);
    opts.parallelWorkers = 1; %parallel CPU matlab toolbox flag; if zero or one no parallelization is used, if a positive integer will attempt to open this number of workers using parallel toolbox
elseif ~isnumeric(opts.parallelWorkers)
    logwrite('Warning: Bad parallelWorkers option. Setting to default (1 = no parallelization).',1);
    opts.parallelWorkers = 1; %parallel CPU matlab toolbox flag; if zero or one no parallelization is used, if a positive integer will attempt to open this number of workers using parallel toolbox
end

if opts.parallelWorkers > 1
    if ~isfield(opts,'closeParallelPool')
        logwrite('No closeParallelPool option found. Setting to default (1 = close parallel pool after finishing)',1);
        opts.closeParallelPool = 1; %close parallel pool after current run? (0 = leave open, 1 = close pool)
    elseif isempty(find([0 1] == opts.closeParallelPool,1))
        logwrite('Warning: Bad closeParallelPool option. Setting to default (1 = close parallel pool after finishing)',1);
        opts.closeParallelPool = 1; %close parallel pool after current run? (0 = leave open, 1 = close pool)
    end
end
    
    
if ~isfield(opts,'savePreProcessed')
    logwrite('No savePreProcessed option found. Setting to default (0 = no)',1);
    opts.savePreProcessed = 0; % save preprocessed data? (0 = no, 1 = yes)
elseif isempty(find([0 1] == opts.savePreProcessed,1))
    logwrite('Warning: Bad savePreProcessed option. Setting to default (0 = no)',1);
    opts.savePreProcessed = 0; % save preprocessed data? (0 = no, 1 = yes)
end

if opts.savePreProcessed == 1
    if ~isfield(opts,'preProcessedSuffix')
        logwrite('No preProcessedSuffix option found. Setting to default (_preprocessed)',1);
        opts.preProcessedSuffix = '_preprocessed';
    elseif ~isstr(opts.preProcessedSuffix)
        logwrite('Bad preProcessedSuffix option. Setting to default (_preprocessed)',1);
        opts.preProcessedSuffix = '_preprocessed';
    end
end

if ~isfield(opts,'saveProcessNetwork')
    logwrite('No saveProcessNetwork option found. Setting to default (1 = yes).',1);
    opts.saveProcessNetwork = 1; % save R structure output from ProcessNetwork computations. 
elseif isempty(find([0 1] == opts.saveProcessNetwork,1))
    logwrite('Warning: Bad saveProcessNetwork option. Setting to default (1 = yes).',1);
    opts.saveProcessNetwork = 1;  % save R structure output from ProcessNetwork computations. 
end

if opts.saveProcessNetwork == 1
    if ~isfield(opts,'outFileProcessNetwork')
        logwrite('No outFileProcessNetwork option found. Setting to default (current date and time).',1);
        opts.outFileProcessNetwork = NaN; % output file will be named with date and time of processing run
    elseif ~ischar(opts.outFileProcessNetwork) && ~isnan(opts.outFileProcessNetwork)
        logwrite('Warning: Bad outFileProcessNetwork option. Setting to default (current date and time).',1);
        opts.outFileProcessNetwork = NaN; % output file will be named with date and time of processing run
    end
end

if opts.saveProcessNetwork == 1 || opts.savePreProcessed == 1

    if ~isfield(opts,'outDirectory')
        logwrite('No saveDirectory option found. Setting to default (current directory).',1);
        opts.outDirectory = '.\'; % save output files to matlab format for subsequent processing
    end
    if ~exist(opts.outDirectory,'dir')
        logwrite(['Warning: The output directory ' opts.outDirectory ' does not exist. Setting outDirectory to current directory.'],1);
        opts.outDirectory = '.\'; % save output files to matlab format for subsequent processing
    end
    if ~strcmp(opts.outDirectory(end),'\')
        opts.outDirectory = [opts.outDirectory '\'];
    end
end