function [H] = multiCouplingSynchronyPlot(R,popts)
%
% Horizontal bar plot showing the zero-lag coupling statistic (black bar) 
% as well as the maximum coupling statistic at any lag (extension to the 
% bar colored according to lag) for couplings between one variable (TO 
% variable) and several others (FROM variables).
% 
% Usage: [H] = multiCouplingSynchronyPlot(R,popts)
%
% ----- Inputs -----
% R = the output structure from a successful run of the ProcessNetwork
%     Software
% popts = a structure of options for plotting. Some are required while
%        others are available as defaults (indicated below). For optional
%        fields, simply don't include the field in the popts structure.
%
% Fields in popts structure (popts.xxx, where xxx is one of the names below)
% NOTE: structure names are case sensitive
% testStatistic = a string indicating a field in the R structure to plot
%       according to lag (this means it must vary according to lag)
% SigThresh (optional) = a string indicating a field in the R structure
%       corresponding with the statistical significance threshold of
%       testStatistic. If absent, no significance threshold will be
%       plotted.
% ToVar = a single integer or string indicating the 'TO' variable in the
%       coupling (i.e. coupling FROM X TO Y, or X > Y). If an integer,
%       ToVar is the index within R.varNames of the chosen variable. If a
%       string, ToVar is the name of the variable matching one of the
%       values within R.varNames.
% FromVars (optional) = a n-element numeric vector or cell array of strings
%       indicating the 'FROM' variables of the desired couplings to plot 
%       (see description of ToVar above). If a cell array, each 
%       string is a variable name within R.varNames. If a numeric vector, 
%       list the indices of the variables within R.varNames. If absent, all
%       variables except the ToVar will be included.
% fi(optional) = an integer indicating the index of the 4th dimension of 
%       the test statistic corresponding to the fi-th processed file. If 
%       absent, the 1st index of the 4th dimension will be used.
% laglim (optional) = a 2-element vector indicating the [min max] lag to
%       evaluate. If absent, the entire R.lagVect will be evaluated.
% claglim (optional) = a 2-element vector indicating the [min max] lag
%       limits for the colorbar (which shows the lag of the max test
%       statistic). If absent, the min and max from the data will be used.
% fignum (optional) = an integer indicating the figure number to plot in.
%       If absent, figure 1 will be used.
% subplot (optional) = a 3-element row vector of integers indicating the
%       [nrows ncolumns index] of the subplot to plot within. This syntax 
%       follows the subplot() function in matlab.
% clearplot (optional) = a binary 0 (no) or 1 (yes) indicating whether to
%       clear the axes before plotting.
% saveFig (optional) = binary 0 (no) or 1 (yes) to save the figure.
% figName (optional) = a string containing the file name to save the
%       figure, without extension. If absent, then a default name will be
%       used including the name of the test statistic, the file index
%       (popts.fi) and the "TO" variable name in the coupling.
% saveFormat (optional) = a string giving the format to save the figure as
%       (without the .). For example: 'png'. See matlab saveas function 
%       for allowable image formats. If absent, the figure will be saved 
%       in the native .fig format.
% outDirectory (optional) = a string giving the full or relative path of
%       the directory to save the figure in. Be sure to include the ending
%       slash (.\myDirectory\). If absent or improper, the figure will be
%       saved in the current directory.
%
% ------ Outputs ------
% fig = the figure handle of the plot
% axis = the axes handle of the plot
% bar0 = the handle(s) of the plotted zero-lag bars
% barMax = the handle(s) of the plotted maximum statistic bars
% sigThresh = the handle(s) of the plotted statistical significance lines
% leg = the handle of the legend
% cb = the handle of the colorbar
%
% --------------------
% Copyright Cove Sturtevant, 2015. 

% Initialize
H.fig = NaN; % figure handle
H.axis = NaN; % axis handle
H.bar0 = NaN; % bar handles for 0-lag statistic
H.barMax = NaN; % bar handles for max statistic
H.sigThresh = NaN; % line handles for statistical significance thresholds
H.leg = NaN; % legend handle

% Check options
if ~isfield(popts,'testStatistic')
    disp('No test statistic indicated. Check options structure.')
    return
elseif ~isfield(R,popts.testStatistic)
    disp([popts.testStatistic ' is not a valid statistic in the results structure.'])
    return
end

if ~isfield(popts,'offset')
    popts.offset = 0;
elseif numel(popts.offset) ~= 1 && (size(popts.offset,1) ~= size(X,1) || size(popts.offset,2) ~= size(popts.offset,2))
    disp('Invalid offset. Check popts.offset')
    return
end

% Pull the test statistic
X = eval(['R.' popts.testStatistic])+popts.offset;
[nVars,~,nLags,nFiles] = size(X);

if ~isfield(popts,'SigThresh')
    XSigThresh = [];
elseif ~isfield(R,popts.SigThresh)
    XSigThresh = [];
else
    XSigThresh = eval(['R.' popts.SigThresh])+popts.offset;    
end

if ~isfield(popts,'ToVar')
    disp('Missing index of ''TO'' variable.')
    return
elseif isnumeric(popts.ToVar) && (min(popts.ToVar) < 1 || max(popts.ToVar) > nVars)
    disp('Invalid ''TO'' variable index.')
    return
elseif iscell(popts.ToVar)
    % Get variable ex
    ci = find(strcmpi(popts.ToVar(1),R.varNames) == 1,1);
    
    if isempty(ci)
        disp('Invalid ''TO'' variable name.')
        return
    end
elseif ischar(popts.ToVar)
    % Get variable ex
    ci = find(strcmpi(popts.ToVar(1,:),R.varNames) == 1,1);
    
    if isempty(ci)
        disp('Invalid ''TO'' variable name.')
        return
    end
else
    ci = popts.ToVar;
end
    
if ~isfield(popts,'FromVars')
    ri = setxor(1:nVars,ci);
elseif isempty(popts.FromVars)
    ri = setxor(1:nVars,ci);    
elseif isnumeric(popts.FromVars) && (min(popts.FromVars) < 1 || max(popts.FromVars) > nVars)
    disp('Invalid ''FROM'' variable indices.')
    return
elseif iscell(popts.FromVars)
    % Get variable indices 
    ri = NaN(length(popts.FromVars),1);
    for i = 1:length(popts.FromVars)
        r = find(strcmpi(popts.FromVars(i),R.varNames) == 1,1);
        if isempty(r)
            disp('At least one ''FROM'' variable name does not match R.varNames'.')
            return
        else
            ri(i) = r;
        end
    end
else
    ri = popts.FromVars;
end

if ~isfield(popts,'fi')
    popts.fi = 1;
elseif popts.fi < 1 || popts.fi > nFiles
    disp('Invalid file index (4th dimension).')
    return
end    

if ~isfield(popts,'laglim')
    popts.laglim = [min(R.lagVect) max(R.lagVect)];
elseif isempty(popts.laglim )
    popts.laglim = [min(R.lagVect) max(R.lagVect)];    
elseif numel(popts.laglim) ~= 2
    popts.laglim = [min(R.lagVect) max(R.lagVect)];
elseif popts.laglim(2) - popts.laglim(1) <= 0
    popts.laglim = [min(R.lagVect) max(R.lagVect)];
end    

if ~isfield(popts,'claglim')
    popts.claglim = [];
elseif isempty(popts.claglim)
    popts.claglim = [];    
elseif numel(popts.claglim) > 2 || popts.claglim(2)-popts.claglim(1) <= 0
    popts.claglim = [];
end    

if ~isfield(popts,'fignum')
    popts.fignum = 1;
elseif isempty(popts.fignum)
    popts.fignum = 1;
elseif ~isnumeric(popts.fignum) || popts.fignum < 0 || isnan(popts.fignum)
    popts.fignum = 1;
end

if ~isfield(popts,'subplot')
    popts.subplot = [];
elseif ~isnumeric(popts.subplot) || numel(popts.subplot) ~= 3 || isnan(sum(popts.subplot))
    disp('Bad subplot option. Plotting as single plot.')
    popts.subplot = [];
end

if ~isfield(popts,'clearplot')
    popts.clearplot = 0;
end

if ~isfield(popts,'saveFig')
    popts.saveFig = 0;
end
if popts.saveFig ~= 0
    if ~isfield(popts,'outDirectory')
        popts.outDirectory = '.\';
    elseif ~exist(popts.outDirectory,'dir')
        disp(['Warning: The output directory ' popts.outDirectory ' does not exist. Setting outDirectory to current directory.'],1);
        popts.outDirectory = '.\';
    end
    
    if ~isfield(popts,'saveFormat')
        popts.saveFormat = 'fig';
    end
end

% Set up figure
figure(popts.fignum); H.fig = gcf;

if ~isempty(popts.subplot)
    subplot(popts.subplot(1),popts.subplot(2),popts.subplot(3))
    H.axis = gca;
    if popts.clearplot ~= 1
        hold on
    else
        delete(H.axis)
        H.axis = subplot(popts.subplot(1),popts.subplot(2),popts.subplot(3));
    end
else
    if popts.clearplot ~= 1
        subplot(1,1,1)
        H.axis = gca;
        hold on
    else
        clf
    end
end

% Restrict data to chosen lag range and couplings
lagi = find(R.lagVect >= popts.laglim(1) & R.lagVect <= popts.laglim(2));
lagVect = R.lagVect(lagi);
X = flipud(reshape(X(ri,ci,lagi,popts.fi(1)),length(ri),length(lagi)));
l0i = find(lagVect == 0,1);

if isempty(l0i)
    disp('ERROR: Lag of 0 not found within chosen lag range. Check options and results.')
    return
end

X0 = X(:,l0i); % zero-lag statistic
[XM,maxi] = max(X,[],2); % max statistic at any lag
lagMax = lagVect(maxi); % lag of max statistic

% Get statistical significance threshold
if ~isempty(XSigThresh)
    if size(XSigThresh,3) ~= nLags
        XsT = flipud(XSigThresh(ri,ci,1,popts.fi(1)));
    else
        XsT = flipud(XSigThresh(ri,ci,maxi,popts.fi(1))); % sig thresh at lag of max statistic
        %XsT = flipud(mean(XSigThresh(ri,ci,:,popts.fi(1)),3));% mean of sig threshes at any lag

    end
else
    XsT = NaN(size(X0));
end

% Assign colors to the Max Statistic according to the lag
cmap = colormap(redblue(100));
if ~isempty(popts.claglim)
    cax = popts.claglim;
elseif numel(lagMax) == 1
    cax = [min([0 min(lagMax) max(lagMax)]) max([0 min(lagMax) max(lagMax)])];
else
    cax = [min(lagMax) max(lagMax)];
end
if cax(2)-cax(1) <= eps
    cax = [min(lagVect) max(lagVect)];
end
cmi = dsearchn(linspace(cax(1),cax(2),100)',lagMax);
colors = cmap(cmi,:);

for i = 1:length(ri)
    % Plot max test statistic colored according to the lag
    H.barMax(i,1) = barh(i,XM(i),0.8,'facecolor',colors(i,:));
    hold on
    
    % Plot test statistic at zero lag
    H.bar0(i,1) = barh(i,X0(i),0.8,'facecolor','k');
    
    % Plot statistical significance threshold
    H.sigThresh(i,1) = plot(XsT(i)*ones(1,2),[i-0.5 i+0.5],'--','color',[0.5 0.5 0.5]);
end
set(gca,'ytick',1:length(ri),'yticklabel',rot90(R.varSymbols(ri),2))
cbh = colorbar;
caxis(cax)
xlabel([popts.testStatistic])
ylabel('Variable Symbol')
set(get(cbh,'ylabel'),'String',['Lag (\tau) of Max ' popts.testStatistic])
set(gca,'ylim',[0 length(ri)+1])
llh = barh(NaN,1,0.8,'facecolor',cmap(50,:));
if ~isempty(XSigThresh)
    H.leg = legend([H.bar0(1);llh;H.sigThresh(1)],[popts.testStatistic '_{X>' R.varSymbols{ci} '}(\tau = 0)'],['Max ' popts.testStatistic '_{X>' R.varSymbols{ci} '}(\tau \neq 0)'],'sigThresh','location','best');
else
    H.leg = legend([H.bar0(1);llh],[popts.testStatistic '_{X>' R.varSymbols{ci} '}(\tau = 0)'],['Max ' popts.testStatistic '_{X>' R.varSymbols{ci} '}(\tau \neq 0)'],'location','best');
end
title([R.varSymbols{ci} ' (file ' num2str(popts.fi) ')'])
set(gca,'box','on')
H.cb = cbh;

if popts.saveFig
    
    if isfield(popts,'figName')
        fname = popts.figName;
    else
        fname = [popts.testStatistic '_' num2str(popts.fi) '_' R.varNames{ci}];
    end
    
    try
        saveas(H.fig,[popts.outDirectory fname],popts.saveFormat);
    catch
        disp(['Problem saving figure: ' [popts.outDirectory fname '.' popts.saveFormat] '. Check options.'])
    end
end

