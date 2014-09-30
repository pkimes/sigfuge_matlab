function SigFugePlotPK(data, paramstruct) 
% SIGFUGEPLOTPK, plot for SigFuge curve view of multiple RNA-seq samples
%                along gene model.
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
%     ipval            whether to calculate SigFuge p-value,
%                       default is 1 (calculate)
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
%last updated: 02/18/2014


%specify default values for parameter structure entries
genename = '';
labels = 0;
ipval = 1;
savestr = 'SigFugePlot';
exonpath = 0;

%use specified parameter values if passed
if nargin > 1;

    if isfield(paramstruct, 'genename');
        genename = paramstruct.genename;
    end;
    
    if isfield(paramstruct, 'labels');
        labels = paramstruct.labels;
    end;

    if isfield(paramstruct, 'ipval');
        ipval = paramstruct.ipval;
    end;

    if isfield(paramstruct, 'exonpath');
        exonpath = paramstruct.exonpath;
    end;

    if isfield(paramstruct, 'savestr');
        savestr = paramstruct.savestr; 
    end;
    
end;


[~, n] = size(data);


%calculate SigFuge labels if not specified
if length(labels) ~= n;
    labels = SigFugeLabelsPK(data);
    disp('!!!!       using SigFuge labels       !!!!');
end;


%calcuate SigFuge p-value
if ipval;
    [~, pvalZ] = SigFugePvalPK(data, labels);
end ;


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
%compute cluster medians
medn1 = median(logdata2(:,labels==1), 2);
medn2 = median(logdata2(:,labels==2), 2);
medn3 = median(logdata2(:,labels==3), 2);

%compute various plotting dimensions
vaxh2 = axisSM(1:size(logdata2, 1));
vaxv2 = axisSM(logdata2);  
mmin2 = min(reshape(logdata2, 1, []));
mmax2 = max(reshape(logdata2, 1, []));
vaxv2(1) = mmin2 - .08*(mmax2 - mmin2);
vaxv2B5  = mmin2 - .06*(mmax2 - mmin2);


%initialize and fill figure/plot
figure(1);
clf;
hold on;
    
plot(logdata2(:,labels==1), '-', 'color', [.4 .4 .4]);
plot(logdata2(:,labels==2), '-', 'color', [1 .3 .3]);
plot(logdata2(:,labels==3), '-', 'color', [.3 .3 1]);

plot(medn1, '-', 'color', [.0 .0 .0], 'linewidth', 3);
plot(medn2, '-', 'color', [.6 .0 .0], 'linewidth', 3);
plot(medn3, '-', 'color', [.0 .0 .6], 'linewidth', 3);

xlabel('exonic nt number, not genomic position');
ylabel('log_{10} (counts+1)');
title([genename ' curves colored by group labels']);

axis([vaxh2(1) vaxh2(2) vaxv2(1) vaxv2(2)]);
ax = axis;


%add exon boundaries if annotations available
if exist('annots', 'var');
    addExBd(annots, 0, vaxv2B5, .03*(mmax2 - mmin2));
end;


%print p-value if provided in paramstruct
if ipval;
  text(ax(1)+(ax(2)-ax(1))*.9, ...
       ax(3)+(ax(4)-ax(3))*1.04, ...
       ['SigFuge p-val = ' num2str(pvalZ)], ...
        'FontSize',10) ;
end ;


%print number of observations per cluster
text(ax(1)+(ax(2)-ax(1))*.9, ...
     ax(3)+(ax(4)-ax(3))*.99, ...
     '#(Red,Blue,Black)/Total:', ...
     'FontSize',10);
text(ax(1)+(ax(2)-ax(1))*.9, ...
     ax(3)+(ax(4)-ax(3))*.96, ...
     ['(' num2str(sum(labels==2)) ',' ...
          num2str(sum(labels==3)) ',' ...
          num2str(sum(labels==1)) ')' ...
      ' / ' num2str(n)], ...
      'FontSize',10);


%add strandedness if annotations available
if exist('annots', 'var');
    text(ax(1)+(ax(2)-ax(1))*.9, ...
         ax(3)+(ax(4)-ax(3))*.92, ...
         ['gene on ' annots.dir ' strand'], ...
         'FontSize',8);
end;


%write to file
orient landscape;
print('-dpsc2', savestr, '-append');



