clear

sitename = 'NZKop';
wavedir = ['C:\Users\csturtevant\Dropbox\KnoxLab\CoveAnalysis\WaveletDecomposed\'];
datafileprefix = [sitename '_ProcessNetworkVariables_raw_withRWSurr'];

testscales = 7:10;

for si = testscales
    
    % Open the detail file
    load([wavedir datafileprefix '_detail_' num2str(si) '.mat'])
    
    % Get Variances
    stdData=std(Data,0,1,'omitnan');
    stdSurr=std(Surrogates,0,1,'omitnan');
    
    Surrogates=Surrogates./repmat(stdSurr,size(Data,1),1,1).*repmat(stdData,size(Data,1),1,1);
    
    save([wavedir datafileprefix '_detail_' num2str(si) '_surrVarMatch.mat'],'Data','Surrogates')
    
    clear Data Surrogates
end