% Imports the retinal sensitivities 
function S = import_RetinalSensitivities(path)

    cd(path)
            
    % common settings
    bands = 'both';
    linLog = 'lin';
    quantaE = 'Q';
    xRes = 1;
    xLims = [380 780];

    S.wavelength = (xLims(1):xRes:xLims(2))';
    
    fieldNames = {'SWS'; 'MWS'; 'LWS'; 'Rhodop'; 'Melanop'};

    % SWS
    k = 1;
    S.(fieldNames{k}).peak_nm = 420; % http://www.ncbi.nlm.nih.gov/pubmed/6140680
    S.(fieldNames{k}).density = 0.30; % http://dx.doi.org/10.1016/S0042-6989(98)00225-9
    S.(fieldNames{k}).data = nomog_Govardovskii2000(S.(fieldNames{k}).peak_nm, bands, linLog, quantaE, S.(fieldNames{k}).density, xRes, xLims);

    % MWS
    k = 2;
    S.(fieldNames{k}).peak_nm = 530; % http://www.ncbi.nlm.nih.gov/pubmed/6140680
    S.(fieldNames{k}).density = 0.38; % http://dx.doi.org/10.1016/S0042-6989(98)00225-9
    S.(fieldNames{k}).data = nomog_Govardovskii2000(S.(fieldNames{k}).peak_nm, bands, linLog, quantaE, S.(fieldNames{k}).density, xRes, xLims);

    % LWS
    k = 3;
    S.(fieldNames{k}).peak_nm = 560; % http://www.ncbi.nlm.nih.gov/pubmed/6140680
    S.(fieldNames{k}).density = 0.38; % http://dx.doi.org/10.1016/S0042-6989(98)00225-9
    S.(fieldNames{k}).data = nomog_Govardovskii2000(S.(fieldNames{k}).peak_nm, bands, linLog, quantaE, S.(fieldNames{k}).density, xRes, xLims);

    % Rhodopsin
    k = 4;
    S.(fieldNames{k}).peak_nm = 495;
    S.(fieldNames{k}).density = 0.40; % http://dx.doi.org/10.1016/0042-6989(95)00114-F
    S.(fieldNames{k}).data = nomog_Govardovskii2000(S.(fieldNames{k}).peak_nm, bands, linLog, quantaE, S.(fieldNames{k}).density, xRes, xLims);

    % Melanopsin
    k = 5;
    S.(fieldNames{k}).peak_nm = 480;
    S.(fieldNames{k}).density = 0.001; % http://dx.doi.org/10.1038/nature07682
    S.(fieldNames{k}).data = nomog_Govardovskii2000(S.(fieldNames{k}).peak_nm, bands, linLog, quantaE, S.(fieldNames{k}).density, xRes, xLims);

        % normalize
        for k = 1 : length(fieldNames)
            % S.(fieldNames{k})
            S.(fieldNames{k}).data = S.(fieldNames{k}).data ./ max(S.(fieldNames{k}).data);
        end

    S.headers = {['SWS (', num2str(S.(fieldNames{1}).peak_nm), ' nm)'];...
                  ['MWS (', num2str(S.(fieldNames{2}).peak_nm), ' nm)'];...
                  ['LWS (', num2str(S.(fieldNames{3}).peak_nm), ' nm)'];...
                  ['Rod (', num2str(S.(fieldNames{4}).peak_nm), ' nm)'];...
                  ['OPN4 (', num2str(S.(fieldNames{5}).peak_nm), ' nm)']};