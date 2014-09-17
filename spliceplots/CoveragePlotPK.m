function CoveragePlotPK(data, paramstruct) 
% COVERAGEPLOTPK, plot for curve view of multiple RNA-seq samples
%                along gene model.
%                   - essentially SigFugePlotPK for general K clusters
%
%
% Inputs:
%   data         - d x n matrix of read counts
%
%   paramstruct  - a Matlab structure of input parameters
%
%     fields           values
%
%     genename         name of locus to get annotations
%     labels           nx1 vector of 1,2,3's corresponding to sample labels
%                       if not specified, will use SigFuge labels,
%                       default = 0
%     exonpath         where to find information on exon annotations,
%                       default is set to 0 and will use default exonpath
%                       for readinRefSeq.m, if no exonpath is to be used,
%                       specify the empty string ''
%     savestr          file name for saving,
%                       default is : 'SigFugePlot'
%
% Output:
%     SFplot     - SigFuge curve plot with necessary information
%
%written by: Patrick Kimes
%last updated: 10/03/2013


%specify default values for parameter structure entries
genename = '';
labels = 0;
medianwidth = 3;
savestr = 'SigFugePlot';
exonpath = 0;
printn = 1;

%use specified parameter values if passed
if nargin > 1;

    if isfield(paramstruct, 'genename');
        genename = paramstruct.genename;
    end;
    
    if isfield(paramstruct, 'labels');
        labels = paramstruct.labels;
    end;

    if isfield(paramstruct, 'medianwidth');
        medianwidth = paramstruct.medianwidth;
    end;
    
    if isfield(paramstruct, 'printn');
        printn = paramstruct.printn;
    end;

    if isfield(paramstruct, 'exonpath');
        exonpath = paramstruct.exonpath;
    end;

    if isfield(paramstruct, 'savestr');
        savestr = paramstruct.savestr; 
    end;
end;


[d, n] = size(data);


%calculate SigFuge labels if not specified
if length(labels) ~= n;
    disp('!!!!       using SigFuge labels       !!!!');
    labels = SigFugeLabelsPK(data);
end;


%obtain gene exon boundaries and direction for plotting                
%load in exon boundaries
if ischar(exonpath) && ...
    ~isempty(genename);

    if isempty(exonpath);
        annots = findinRefSeq({genename});
    else
        annots = findinRefSeq({genename}, exonpath);
    end;

end;

%flip data orientation if 
if exist('annots', 'var');
    if annots.dir == '-';
        data = flipud(data);
    end;
end;

logdata2 = log10(data+1);

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~iscell(labels);
    labels = textscan(sprintf('%i\n',labels'),'%s');
    labels = labels{1};
end;

tags = unique(labels);
K = length(tags);


medn = zeros(d, K);
%compute cluster medians
for k = 1:K;
    medn(:, k) = median(logdata2(:, logical(strcmp(tags{k}, labels))), 2);
end;

%compute various plotting dimensions
vaxh2 = axisSM(1:size(logdata2, 1));
vaxv2 = axisSM(logdata2);  
mmin2 = min(reshape(logdata2, 1, []));
mmax2 = max(reshape(logdata2, 1, []));
vaxv2(1) = mmin2 - .08*(mmax2 - mmin2);
vaxv2B5  = mmin2 - .06*(mmax2 - mmin2);


cmap = colormap('HSV');
colors = cmap(floor(linspace(1, size(cmap,1), K+1)), :);
colors = colors(1:K, :);

%initialize and fill figure/plot
figure(1);
clf;
hold on;

for k = 1:K;
    plot(logdata2(:,logical(strcmp(tags{k}, labels))), '-', ...
         'color', colors(k,:));
end;

if medianwidth > 0;
    for k = 1:K;
        plot(medn(:,k), '-', ...
             'color', max(colors(k,:)*.75,0), ...
             'linewidth', medianwidth);
    end;
end;

xlabel('exonic nt number, not genomic position');
ylabel('log_{10} (counts+1)');
title([genename ' curves colored by group labels']);

axis([vaxh2(1) vaxh2(2) vaxv2(1) vaxv2(2)]);
ax = axis;


%add exon boundaries if annotations available
if exist('annots', 'var');
    addExBd(annots, 0, vaxv2B5, .03*(mmax2 - mmin2));
end;

%calculate counts
table = tabulate(labels);
printlab = cellfun(@(x) [x ', '], table(:,1), ...
                   'uniformoutput', 0);
printlab = [printlab{:}]; printlab = printlab(1:(end-2));

printcnt =  cellfun(@(x) [num2str(x) ', '], table(:,2), ...
                   'uniformoutput', 0);
printcnt = [printcnt{:}]; printcnt = printcnt(1:(end-2));

%print number of observations per cluster
if printn;
    text(ax(1)+(ax(2)-ax(1))*.6, ...
         ax(3)+(ax(4)-ax(3))*.99, ...
         ['#(' printlab ')/Total:'], ...
         'FontSize',10);
    text(ax(1)+(ax(2)-ax(1))*.6, ...
         ax(3)+(ax(4)-ax(3))*.96, ...
         ['(' printcnt ')' ...
          ' / ' num2str(n)], ...
          'FontSize',10);
    %add strandedness if annotations available
    if exist('annots', 'var');
        text(ax(1)+(ax(2)-ax(1))*.9, ...
             ax(3)+(ax(4)-ax(3))*.92, ...
             ['gene on ' annots.dir ' strand'], ...
             'FontSize',8);
    end;
else
    legend(tags);
end;



%write to file
orient landscape;
print('-dpsc2', savestr, '-append');



