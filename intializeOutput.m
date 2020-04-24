function [R] = intializeOutput(nDataFiles,nVars,opts)
nBins = opts.nBins'; % Make column vector
nLags = length(opts.lagVect);

% Determine size of 3rd dimension for statistical tests depending on 
% whether we test each lag
if opts.SurrogateMethod > 0
    if opts.SurrogateTestEachLag
        nSLags = nLags;
    else
        nSLags = 1;
    end
else
    nSLags = 1;
end

% INITIALIZE PREPROCESSING QUANTITIES
R.nRawData               = zeros(nDataFiles,1);
R.nVars                  = zeros(nDataFiles,1);
R.varNames               = opts.varNames;
R.varSymbols             = opts.varSymbols;
R.varUnits               = opts.varUnits;
if length(nBins) < nVars
    R.nBinVect           = ones(nVars,1)*nBins;
else
    R.nBinVect           = nBins;
end
R.nClassified            = NaN(nDataFiles,1);
R.binEdgesLocal          = NaN(nVars,max(nBins),nDataFiles);
R.minEdgeLocal           = NaN(nVars,nDataFiles);
R.maxEdgeLocal           = NaN(nVars,nDataFiles);
R.minSurrEdgeLocal       = NaN(nVars,nDataFiles);
R.maxSurrEdgeLocal       = NaN(nVars,nDataFiles);
R.LocalVarAvg            = NaN(nVars,nDataFiles);
R.LocalVarCnt            = NaN(nVars,nDataFiles);
R.binEdgesGlobal         = NaN(nVars,max(nBins));
R.minEdgeGlobal          = NaN(nVars,1);
R.maxEdgeGlobal          = NaN(nVars,1);
R.binSurrEdgesGlobal     = NaN(nVars,max(nBins));
R.minSurrEdgeGlobal      = NaN(nVars,1);
R.maxSurrEdgeGlobal      = NaN(nVars,1);
R.GlobalVarAvg           = NaN(nVars,1);

% INITIALIZE OUTPUT QUANTITIES FROM THE ENTROPYFUNCTION
R.lagVect                = opts.lagVect;
R.HXt                    = NaN(nVars,nVars,nLags,nDataFiles);
R.HYw                    = NaN(nVars,nVars,nLags,nDataFiles);
R.HYf                    = NaN(nVars,nVars,nLags,nDataFiles);
R.HXtYw                  = NaN(nVars,nVars,nLags,nDataFiles);
R.HXtYf                  = NaN(nVars,nVars,nLags,nDataFiles);
R.HYwYf                  = NaN(nVars,nVars,nLags,nDataFiles);
R.HXtYwYf                = NaN(nVars,nVars,nLags,nDataFiles);
R.SigThreshT             = NaN(nVars,nVars,nSLags,nDataFiles);
R.SigThreshI             = NaN(nVars,nVars,nSLags,nDataFiles);
R.SigThreshL             = NaN(nVars,nVars,nSLags,nDataFiles);
R.meanShuffT             = NaN(nVars,nVars,nSLags,nDataFiles);
R.sigmaShuffT            = NaN(nVars,nVars,nSLags,nDataFiles);
R.meanShuffI             = NaN(nVars,nVars,nSLags,nDataFiles);
R.sigmaShuffI            = NaN(nVars,nVars,nSLags,nDataFiles);
R.meanShuffL             = NaN(nVars,nVars,nSLags,nDataFiles);
R.sigmaShuffL            = NaN(nVars,nVars,nSLags,nDataFiles);
R.nCounts                = NaN(nVars,nVars,nLags,nDataFiles);
R.I                      = NaN(nVars,nVars,nLags,nDataFiles);
R.T                      = NaN(nVars,nVars,nLags,nDataFiles);
R.L                      = NaN(nVars,nVars,nLags,nDataFiles); % ADDED  linear redundancy
R.IR                     = NaN(nVars,nVars,nLags,nDataFiles);
R.TR                     = NaN(nVars,nVars,nLags,nDataFiles);
R.SigThreshTR            = NaN(nVars,nVars,nSLags,nDataFiles);
R.SigThreshIR            = NaN(nVars,nVars,nSLags,nDataFiles);
R.SigThreshIRnormByShuff = NaN(nVars,nVars,nSLags,nDataFiles); % ADDED significance threshold for IR normalized by meanShuffIR
R.meanShuffTR            = NaN(nVars,nVars,nSLags,nDataFiles);
R.sigmaShuffTR           = NaN(nVars,nVars,nSLags,nDataFiles);
R.meanShuffIR            = NaN(nVars,nVars,nSLags,nDataFiles);
R.sigmaShuffIR           = NaN(nVars,nVars,nSLags,nDataFiles);
R.IRnormByShuff          = NaN(nVars,nVars,nSLags,nDataFiles); % ADDED IR normalized (divided by) meanShuffIR
R.sigmaShuffIRnormByShuff= NaN(nVars,nVars,nSLags,nDataFiles); % ADDED signmaShuff IR normalized (divided by) meanShuffIR
R.IRsubtractShuff        = NaN(nVars,nVars,nSLags,nDataFiles); % ADDED IR normalized (subtracting out) meanShuffIR
R.IRSigThreshIRsubtractShuff= NaN(nVars,nVars,nSLags,nDataFiles); % ADDED sig thresh for IRsubtractShuff
R.Tplus                  = NaN(nVars,nLags,nDataFiles);
R.Tminus                 = NaN(nVars,nLags,nDataFiles);
R.Tnet                   = NaN(nVars,nLags,nDataFiles);
R.TnetBinary             = NaN(nVars,nVars,nLags,nDataFiles);
R.InormByDist            = NaN(nVars,nVars,nLags,nDataFiles);
R.TnormByDist            = NaN(nVars,nVars,nLags,nDataFiles);
R.SigThreshInormByDist   = NaN(nVars,nVars,nSLags,nDataFiles);
R.SigThreshTnormByDist   = NaN(nVars,nVars,nSLags,nDataFiles);
R.Ic                     = NaN(nVars,nVars,nLags,nDataFiles);
R.Tc                     = NaN(nVars,nVars,nLags,nDataFiles);
R.TvsIzero               = NaN(nVars,nVars,nLags,nDataFiles);
R.SigThreshTvsIzero      = NaN(nVars,nVars,nSLags,nDataFiles);
R.IvsIzero               = NaN(nVars,nVars,nLags,nDataFiles);
R.SigThreshIvsIzero      = NaN(nVars,nVars,nSLags,nDataFiles);
R.Abinary                = NaN(nVars,nVars,nLags,nDataFiles);
R.Awtd                   = NaN(nVars,nVars,nLags,nDataFiles);
R.AwtdCut                = NaN(nVars,nVars,nLags,nDataFiles);
R.charLagFirstPeak       = NaN(nVars,nVars,nDataFiles);
R.TcharLagFirstPeak      = NaN(nVars,nVars,nDataFiles);
R.charLagMaxPeak         = NaN(nVars,nVars,nDataFiles);
R.TcharLagMaxPeak        = NaN(nVars,nVars,nDataFiles);
R.TvsIzerocharLagMaxPeak = NaN(nVars,nVars,nDataFiles);
R.nSigLags               = NaN(nVars,nVars,nDataFiles);
R.FirstSigLag            = NaN(nVars,nVars,nDataFiles);
R.LastSigLag             = NaN(nVars,nVars,nDataFiles);
R.HXtNormByDist          = NaN(nVars,nVars,nLags,nDataFiles);
