function annotations = findinRefSeq(genenames, exonpath)
%
%
% Inputs:
%   genenames  - cell array of gene names
%   exonpath   - path of exon_boundaries.txt file, optional
%                 default is to use UCSC hg19 annotations from 
%                 TCGA GAF v2.1
%   
% Outputs:
%   annotations - structure with entries
%                   genename, 
%                   dir,
%                   (exon) start, 
%                   (exon) stop
%
%
%written by: Patrick Kimes
%last updated: 10/02/2013


if nargin > 1;
    RefSeq = readinRefSeq(exonpath);
else
    RefSeq = readinRefSeq();
end;

failed = false;

if iscell(genenames);
    
    if all(cellfun(@ischar, genenames));
        findidx = @(x) find(strcmp(x, {RefSeq.gene}));
        idx = cellfun(findidx, genenames);
    else
        failed = true;
    end;
    
elseif ischar(genenames);
    
    idx = find(strcmp(genenames, {RefSeq.gene}));

else
    
    failed = true;

end;

%what to do if failed anywhere above
if failed; %if input wasn't a cell array or string
    
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
    disp('!!! input must be a gene name or cell array !!');
    disp('!!! of gene names, breaking without output  !!');
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
    return;   
    
end;

    %pull out subset of annotations
    annotations = RefSeq(idx);
    
    %convert string of exon start/ends delimited by commas to vector
    %doing this using an anonymous function
    getit = @(x) cellfun(@str2num, regexp(x, ',', 'split'));
    for i = 1:length(genenames);
        annotations(i).start = getit(annotations(i).start);
    end;

    for i = 1:length(genenames);
        annotations(i).stop = getit(annotations(i).stop);
    end;

    %remove RefSeq from memory/workspace
    clear RefSeq;
    
end

