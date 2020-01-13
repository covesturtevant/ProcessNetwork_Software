% Setup file to run the ProcessNetwork software, and a couple plotting
% options

%% Create the options structure
clear
close all
clc
clear global processLog
global processLog

% File specification and variable designation
opts.files = {'ChaosInputData\data_ChaosWithResolutions.mat'}; % cell column of strings listing files to process, including path
opts.varNames = {'A','B'}; % cell row of strings listing the variable names (no latex symbols here) (default = V1,V2,V3,etc.)
opts.varSymbols = {}; % cell row of strings listing the variable symbols for plotting (latex okay here) (default = same as varNames)
opts.varUnits = {};  % (optional) cell row of strings listing the variable symbols for plotting (latex okay here)

% Preprocessing options
opts.NoDataCode = NaN; % Numerical value representing no data (default = NaN)
opts.trimTheData = 1; % Remove entire rows of data with at least 1 missing value? 0 = no (default), 1 = yes 
opts.transformation = 0; % 0 = apply no transformation (default), 1 = apply anomaly filter, 2 = apply maximal overlap wavelet transform (requires no data gaps)
    opts.anomalyPeriodInData = 24; % Anomaly filter only: set to the length in time steps of the period of the data (default = 48)
    opts.anomalyMovingAveragePeriodNumber = 5; % Anomaly filter only: the moving window used for anomaly generation, in # of periods (default = 5)
    opts.waveN = 1; % Wavelet transform only: vector of dyadic (2^n) scales to decompose (default = 1)
    opts.waveName = 'la8'; % Wavelet transform only: string indicating mother wavelet to use (default = 'la8', options are 'haar','d4','la8','la16')
    opts.waveDorS = 1; % Wavelet transform only: 1 (default) = output detail added at this scale (will be summed over scales in waveN), 2 = output smooth approximation at the scale in waveN (waveN must be single-valued)
opts.binType = 1; % 0 = don't do classification (or data are already classified), 1 = classify using locally bounded bins (default), 2 = classify using globally bounded bins
    opts.binPctlRange = [0 100]; % [min max] percentile range to determine total bin range for each variable. Default is entire data range ([0 100]) 
    opts.nBins = 11; % Number of bins to classify each signal (default = 11). This can be a single value to be applied to all variables, or a row vector of integers correspoding with each column in the input data
opts.SurrogateMethod = 2; % 0 = Do not do statistical testing, regardless of whether input files contain Surrogates. 1 = Use the Surrogates contained in the loaded files. 2 = create and test new surrogates via random shuffle of input data (default); 3 = create and test new surrogates via IAAFT method (requires no data gaps); 4 = create and test new surrogates using random walk method (retains original gaps in data)
    opts.nTests = 100; % Number of surrogates to create and/or test (default = 100)

% Entropy computation options
opts.doEntropy = 1; % run entropy calculations? 0 = no (default), 1 = yes
    opts.lagVect = 0:100; % lags (in units of time steps) to evaluate (default = 0:10). Note: 0 will always be included, whether it is entered here or not.
    opts.SurrogateTestEachLag = 0; % Test surrogates at each lag? 0 = no, test last lag only (default), 1 = yes
    opts.oneTailZ = 1.66; % one-tail z-score for 95% significance given number of tests, 1.66 for 100, 1.68 for 50, 1.71 for 25
opts.parallelWorkers = 2; % parallel CPU matlab toolbox flag; if 0 or 1 no parallelization is used (default), if a positive integer will attempt to open this number of workers using parallel toolbox
    opts.closeParallelPool = 1; % close parallel pool after finishing? 0 = no, 1 = yes (default)

% Saving options
opts.savePreProcessed = 0; % save preprocessed data? 0 = no (default), 1 = yes
    opts.preProcessedSuffix = '_preprocessed'; % string indicating the suffix to add onto the end of each file name when saving preprocessed data
opts.saveProcessNetwork = 1; % save R structure output and process log from ProcessNetwork computations? 0 = no, 1 = yes (default) 
    opts.outFileProcessNetwork = 'Results_ChaosWithResolutions'; % file name of saved ProcessNetwork output. Default file will be named with date and time of processing run
opts.outDirectory = 'ChaosInputData\'; % directory to save output. (default is current directory)


%% Run the script
tic

[R,opts] = ProcessNetwork(opts);

toc

%% Do coupling lag plot
clear popts

popts.testStatistic = 'T'; % String indicating the name of the test statistic within the R. output structure to plot 
    popts.SigThresh = 'SigThreshT'; % (optional) string inficating the name of the statistical significance threshold for the test statistic
popts.vars = [2 1]; % cell array of strings (matching R.varNames) or 2-element vector of indices within R.varNames indicating the variables in the coupling 1st is FROM variable, 2nd is TO variable (default is [1 2])
    popts.fi = 1; % file index (4th dimension of output) (default is 1)
popts.laglim = [0 100]; % [min max] of lags to show (default is entire range)
popts.fignum = 1; % figure number to plot within (default is 1)
    popts.subplot = [1 1 1]; % Following the matlab format [nrows ncols n], indicate the subplot to plot within (default is [1 1 1]
    popts.clearplot = 0; % 1 = clear the current plot before plotting, 0 = plot over whatever is already in the plot (default)
% popts.plotProperties = {'color',[0 0 0],'linewidth',2}; % (optional) following the Name, Value linespec syntax, enclose desired linespec properties for the testStatistic in curly brackets
% popts.sigThreshPlotProperties = {'color',[0 0 0],'linewidth',1};  % (optional) following the Name, Value linespec syntax, enclose desired linespec properties for sigThresh in curly brackets
popts.saveFig = 1; % 1 = save the figure, 2 = don't save (default)
    popts.figName = 'couplingPlot'; % file name for saving figure. Default includes the testStatistic and variables in the coupling
    popts.saveFormat = 'png'; % format to save the figure in. See MATLAB help for saveas function for allowable formats
%     popts.outDirectory = opts.outDirectory; % relative or absolute directory to save figure in (default is current directory)

 % Plot it!
[H] = couplingLagPlot(R,popts);

%% Do multi-coupling synchrony plot
clear popts

popts.testStatistic = 'TR'; % String indicating the name of the test statistic within the R. output structure to plot 
    popts.SigThresh = 'SigThreshTR'; % (optional) string inficating the name of the statistical significance threshold for the test statistic
popts.ToVar = 2; % The 'TO' variable in the set of couplings to plot. Can be a string or index matching one of the names in R.varNames
    popts.FromVars = []; % The 'FROM' variables in the set of couplings to plot. Can be a cell array of strings or vector of indices matching names in R.varNames
    popts.fi = 1; % file index (4th dimension of output) (default is 1)
popts.laglim = [0 100]; % [min max] of lags to evaluate (default is entire range)
    popts.claglim = []; % [min max] range of lags to form the color range in plotting the lag of max testStatistic (default is entire range observed in the data)
popts.fignum = 2; % figure number to plot within (default is 1)
    popts.subplot = [1 1 1]; % Following the matlab format [nrows ncols n], indicate the subplot to plot within (default is [1 1 1]
    popts.clearplot = 1; % 1 = clear the current plot before plotting, 0 = plot over whatever is already in the plot (default)
popts.saveFig = 1; % 1 = save the figure, 2 = don't save (default)
    popts.figName = 'multicouplingPlot'; % file name for saving figure. Default includes the testStatistic and variables in the coupling
    popts.saveFormat = 'png'; % format to save the figure in. See MATLAB help for saveas function for allowable formats
    popts.outDirectory = opts.outDirectory; % relative or absolute directory to save figure in (default is current directory)

 % Plot it!
[H] = multiCouplingSynchronyPlot(R,popts);

