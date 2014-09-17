function [labels] = SigFugeLabelsPK(data, labels, filter)
%SigFugeLabelsPK, produce nx1 vector of labels for dxn input according to
%                 SigFuge clustering rules.
%
% Inputs:
%   data     - d x n matrix of input values (not yet transformed)
%   labels   - n x 1 vector of 1,2 cluster labels
%   filter   - boolean value whether to filter out low expression samples
%              from texting procedure
%
% Outputs:
%   labels   - n x 1 vector of cluster labels 1, 2, 3
%
%
% written by: Patrick Kimes
% last updated: 02/08/2014

%determine size of dataset
[d, n] = size(data);


%determine which procedures need to be carried out on the input data
% based on the input arguments
if (nargin > 1);
    %check validity of labels
    if (length(labels) == n) && ...
       all(unique(labels) == [1 2]');
        if (nargin < 3) || (filter == false);
            disp('labels specified and filtering not desired.');
            disp('just returning labels orginally passed to function');
            disp('(but as labels+1).');
            labels = labels+1;
            return;
        end;
    else 
        disp('labels must be nsamples x 1 vector of 1s and 2s.');
        disp('using usual SigFuge labeling scheme.');
        filter = false;
    end;
else
    
    filter = false;
end;
    

%if function hasn't 'return'ed yet, flag 'low expression' samples
medcov = zeros(n, 1);
percov = zeros(n, 1);
for i = 1:n;
  temp = data(:, i);
  tempflag = (temp == 0);
  percov(i) = 1 - sum(tempflag)/d;
  medcov(i) = median(temp(~tempflag));
end;
medcov(isnan(medcov)) = 0;
flag = (medcov < 5) | (percov < .10);


%if just filtering on input labels is desired,
% we should skip normalization and clustering steps
if filter;
    labels = labels+1;
    labels(flag) = 1;

else 
    %need to see if all samples are above the threshold
    if sum(~flag) < 10;
        vclass = ones(sum(~flag), 1);
    else
        %normalize and log-scale data
        coverage = sum(data, 1);
        datanorm = data(:, ~flag) ./ ...
                    vec2matSM(coverage(~flag)+1, d) * median(coverage(~flag));
        logdata = log10(datanorm+1);

        %apply 2-means clustering on norm'd data
        paramstruct = struct('nrep', 100, ...
                             'iscreenwrite', 0);
        vclass = SigClust2meanRepSM(logdata, paramstruct);

        % re-order labels by first entry as in R/SigFuge
        if vclass(1) == 2;
            vclass = 3-vclass;
        end;
    end;
    
    %combine norm'd data labels and low exp flags
    labels = ones(n, 1);
    labels(~flag) = vclass+1;
end;


end

