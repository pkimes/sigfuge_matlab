function RefSeq = readinRefSeq(exonpath)
%
%script to pull in RefSeq annotations into a single structure with
%the following fields:
%   genename, 
%   dir, 
%   (exon) start, 
%   (exon) stop
%
%currently using UCSC hg19 annotations from TCGA GAF v2.1 as default
%
%written by: Patrick Kimes
%last updated: 10/02/2013


%pull in gene annotations from local file
if nargin < 1;
    exonpath = ['/Users/pkimes/Dropbox/UNC/Statistics/Research/' ...
            'Projects/NextGen_SigFuge/MyCode_LBG/datafiles/' ...
            'exon_boundaries.txt'];
end;

fid = fopen(exonpath);
pat = '%s %s %s %s %s';
output = textscan(fid, pat, ...
                  'Delimiter', '\t', ...
                  'Headerlines', 0, ...
                  'bufsize', 100000);
fclose(fid);

%place exon annotations in nicer gene annotation structure
output = [output{:}];
RefSeq = cell2struct(output, ...
                {'gene', 'dir', 'chr', 'start', 'stop'}, ...
                2);