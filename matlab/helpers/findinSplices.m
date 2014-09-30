function splices = findinSplices(genenames, splicepath, n)
%
%
% Inputs:
%   genenames  - cell array of gene names
%   splicepath - path to splice matrix created by SplicePull_J.pl
%   n          - sample size
%   
% Outputs:
%   annotations - structure with entries
%                   chr,            (string)
%                   dir,            (string)
%                   gene,           (string)
%                   (splice) start, (int)
%                   (splice) stop,  (int)
%                   vals            (double vector)
%
%
%written by: Patrick Kimes
%last updated: 10/13/2013


if nargin < 3;
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
    disp('!!!!!  must specify splicepath, genenames  !!!!!');
    disp('!!!!!  and sample size to work,            !!!!!');
    disp('!!!!!  exiting without output              !!!!!');
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
    return;
end;

%read in splice information
spliceStruct = readinSplices(splicepath, n);


failed = false;

if iscell(genenames);
    
    if all(cellfun(@ischar, genenames));
        findidx = @(x) find(strcmp(x, {spliceStruct.gene}));
        idx = cellfun(findidx, genenames);
    else
        failed = true;
    end;
    
elseif ischar(genenames);
    
    idx = find(strcmp(genenames, {spliceStruct.gene}));

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
    splices = spliceStruct(idx);

    %remove spliceStruct from memory/workspace
    clear spliceStruct;
    
end

