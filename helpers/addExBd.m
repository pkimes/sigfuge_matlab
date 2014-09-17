function addExBd(annotation, xshift, y, height)
%
%script to add exon boundaries based on annotation input, currently
%set to work with annotations drawn from TCGA GAF v2.1.
%
%input:
%   annotation - structure of size 1 (single gene) with entries:
%                   genename, 
%                   dir, 
%                   (exon) start, 
%                   (exon) stop
%   xshift     - horizontal shift (starting pt of exons)
%   y          - bottom of where annotations should be placed
%   height     - how tall the annotation bar should be
%
%written by: Patrick Kimes
%last updated: 09/12/2013

%colors to use for alternating exon boundaries
mcolor = [.1 .3 .7; 0 .5 .6];


exonLens = annotation.stop-annotation.start;

%need to flip order of exons if on negative strand
if annotation.dir == '-';
    exonLens = fliplr(exonLens);
end;


%define exon boundaries
exbd = [0 cumsum(exonLens)] ...
        + xshift;

%plotting exon boundaries for the gene
for i = 2:length(exbd);
    rectangle('Position', [exbd(i-1), y, ...
               exbd(i)-exbd(i-1), height], ...
              'FaceColor', mcolor(1+mod(i,2),:));
end;

end

