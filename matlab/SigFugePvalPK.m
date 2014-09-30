function [pvalQ, pvalZ] = SigFugePvalPK(data, labels, paramstruct)
%SigFugeLabelsPK, produce nx1 vector of labels for dxn input according to
%                 SigFuge clustering rules.
%
% Inputs:
%   data     - d x n matrix of input values (not yet transformed)
%   labels   - n x 1 vector of cluster labels 1,2,3
%   filter
%   paramstruct
%            - parameter structure to pass to SigClust if default SigFuge
%               parameter values are not appropriate.
%               Current default parameters are:
%                  'iCovEst'        2
%                  'iBGSDdiagplot'  0
%                  'nsim'           100
%                  'iCovEdiagplot'  0
%                  'ipValplot'      0
%                  'iscreenwrite'   0
%               look at SigClustSM.m for further details on possible 
%               parameter values 
%
% Outputs:
%   pvalZ    - SigFuge/SigClust Gaussian approximate p-value
%   pvalQ    - SigFuge/SigClust empirical p-value
%
%
% written by: Patrick Kimes
% last updated: 02/08/2014


%if paramstruct is not passed by user, we will
%use our default parameter structure
if nargin < 3;    
    paramstruct = struct('iCovEst', 2, ...
                         'iBGSDdiagplot', 0, ...
                         'nsim', 100, ...
                         'iCovEdiagplot', 0, ...
                         'ipValplot', 0, ...
                         'iscreenwrite', 0, ...
                         'vclass', labels(labels~=1)-1);
                 %note above that vclass must be all 1s and 2s
end;


%only compute SF p-value if enough samples are passed
if sum(labels~=1) >= 5 && ...
    sum(labels==2) > 0 && ...
    sum(labels==3) > 0;
    
    %compute normalized data from passed labels
    d = size(data, 1);
    coverage = sum(data, 1);
    datanorm = data(:, labels~=1) ./ ...
                vec2matSM(coverage(labels~=1)+1, d) ...
                * median(coverage(labels~=1));
    logdata = log10(datanorm+1);

    %actually compute p-values
    [pvalQ, pvalZ] = SigClustSM(logdata, paramstruct);

else
    
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
    disp('!!! not enough samples in labels 2,3  !!!');
    disp('!!! to calculate SigFuge p-values     !!!');
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
    pvalQ = -1; 
    pvalZ = -1;
    
end;

end

