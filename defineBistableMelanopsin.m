function [alphaR_rel, alphaM_rel] = defineBistableMelanopsin(Rpeak,Mpeak,path)

    % create nomograms
    cd(path.nomogram)
    alphaR_rel = nomog_Govardovskii2000(Rpeak);
    alphaM_rel = nomog_Govardovskii2000(Mpeak);            
    cd(path.mainCode)

    % accessory variables, for rate constant calculation
    alpha_multip = 1;
    sigma = 1;
    nm_spacing = 10^-9; 

    