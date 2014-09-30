function [] = WGsigfuge(genelist, coveragepath, savestr, labels)
%WGsigfuge wrapper function for applying SigFuge to a collection of samples
%and genes
%   This function was written to make applying SigFuge to a large
%   collection of samples and genes, e.g. in a Whole Genome (WG) analysis,
%   easier when labels are provided for the samples. 
%   This function can also be used for calculating the p-values for a
%   single gene locus, e.g.
%       if we have a file: /datapath/coverages/CDKN2A_coverage.txt
%
%       >> WGsigfuge('CDKN2A', '/datapath/coverages/');
%
%
% inputs:
%   genelist        - ngenes x 1 column vector of string gene names
%   coveragepath    - either:
%                     1. string path for coverage files.
%                       Individual files are assumed to be of the 
%                       form: [coveragepath genelist{g} '_coverage.txt']
%                       with size (nsamples x d)
%                     2. nsamples x 1 column vector of string bam paths
%                       coverage is then pulled from these nsamples files
%                       for analysis - NOT YET IMPLEMENTED
%   savestr         - string name for output file to be saved. If not
%                      specified, default is 'WGsigfuge_out'.
%   labels          - optional input column vector of 1,2 labels indicating
%                      a priori class labels of interest, e.g. normal (1)
%                      and tumor (2) samples. Vector must be nsamples x 1
%                      or function will terminate without running analysis.
%                      If unspecified, the default SigFuge procedure based
%                      on filtering low samples, normaliation and applying
%                      k-means will be used to determine cluster labels.
%                      
% output:
%   savestr.mat     - a Matlab datafile
%       'pvalQ'        - ngenes x 1 vector of SigFuge Gaussian p-values
%       'pvalZ'        - ngenes x 1 vector of SigFuge empirical p-values
%       'genelist'     - genelist input for the analysis
%       'SFlabels'     - ngenes x nsamples matrix of labels
% 
%
% dependencies:
%   SigFugeLabelsPK.m
%   SigFugePvalPK.m
%
%
% written by: Patrick Kimes
% last updated: 02/08/2014


%if no class labels are specified for the analysis,
% use SigFuge normalization and k-means procedure to 
% obtain class labels
if nargin < 3;
    savestr = 'WGsigfuge_out';
    labels = [];
elseif nargin < 4;
    labels = [];
end;

ngenes = length(genelist);
nsamples = size(coveragepath, 1);
if nsamples == 1;
    disp('will use coverage files already generated at specified path');
    easyfiles = true;
else 
    disp('will pull coverage from specified bam files');
    easyfiles = false;
end;

pvalQ = zeros(ngenes, 1);
pvalZ = zeros(ngenes, 1);

disp(['running SigFuge analysis on: ' savestr]);
for g = 1:ngenes;

    if easyfiles;
        data = textread([coveragepath genelist{g} '_coverage.txt']); %#ok
        data = data';
        nsamples = size(data, 2);
    else 
        disp('unfortunately BAM option not yet implemented...');
        % still need to implement Matlab functions to read coverage/depth
        % from BAM files
    end;
    
    if g==1;
        SFlabels = zeros(nsamples, ngenes);
    end;
    
    %calculate SigFuge labels
    SFlabels(:, g) = SigFugeLabelsPK(data, labels, true);
    %determine whether SigFugeLabelsPK output includes more 
    % than a single non-background cluster
    if max(SFlabels(:, g)) < 3;
        disp('skipping gene for not having enough samples.');
        pvalQ(g) = -1;
        pvalZ(g) = -1;
    else
        [pvalQ(g), pvalZ(g)] = SigFugePvalPK(data, SFlabels(:, g)');
    end;
%     %other option would be to run this without even filtering
%     SFlabels(:, g) = labels+1; %equivalent to below
%     SFlabels(:, g) = SigFugeLabelsPK(data, labels, false);
%     [pvalQ(g), pvalZ(g)] = SigFugePvalPK(data, SFlabels(:, g));

    disp(['done with gene number ' num2str(g) '/' num2str(ngenes)]) ;

end;

save([savestr '.mat'], ...
     'pvalQ', 'pvalZ', 'genelist', 'SFlabels');


end

