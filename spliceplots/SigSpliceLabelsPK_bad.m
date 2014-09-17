function [labels] = SigSpliceLabelsPK(depth,splice,percentSplice)
%SigSpliceLabelsPK, produce nx1 vector of labels for dxn input according to
%                 SigFuge clustering rules.
%                  Function differs from SigFugeLabelsPK by applying the
%                 clustering to depth + splice information.
%
% Inputs:
%   data     - d x n matrix of input values (not yet transformed)
%
% Outputs:
%   labels   - n x 1 vector of cluster labels 1,2,3

[d,n] = size(depth) ;

medcov = zeros(n,1) ;
percov = zeros(n,1) ;
for i = 1:n ;
  temp = depth(:,i) ;
  tempflag = (temp == 0) ;
  percov(i) = 1 - sum(tempflag)/d ;
  medcov(i) = median(temp(~tempflag)) ;
end ;
medcov(isnan(medcov)) = 0 ;
flag = ((medcov < 5) | (percov < .10)) ;

dsplice = size(splice,1) ;

  % calculate normalization factor for equal sample coverage
coverage = sum(depth,1) ;
vscale = (coverage(~flag)+1) / median(coverage(~flag)) ;
%       vscale = vec2matSM(coverage(~flag)+1,d) * median(coverage(~flag)) ;

  % coverage normalize and log transform
depthnorm = depth(:,~flag) ./ vec2matSM(vscale,d) ;
splicenorm = splice(:,~flag) ./ vec2matSM(vscale,dsplice) ;
logdepthn = log10(depthnorm+1) ;
logsplicen = log10(splicenorm+1) ;

  % scale to equal sum of sq. for splices and read depth
TSS1 = sum(n*std(logdepthn,1,2).^2) ;
TSS2 = sum(n*std(logsplicen,1,2).^2) ;

nlogdepth = logdepthn / sqrt(TSS1) ;
nlogsplice = logsplicen / sqrt(TSS2) ;

percentSplice = percentSplice / 100 ;

  % combine data matrices and apply SigClust
if TSS2 == 0 ;
  fullmat = nlogdepth ;
else
  fullmat = [sqrt(1-percentSplice)*nlogdepth ; ...
             sqrt(percentSplice)*nlogsplice] ;
end ;  

paramstruct = struct('nrep',100, ...
                     'iscreenwrite',0) ;
vclass = SigClust2meanRepSM(fullmat,paramstruct) ; 
if sum(vclass==1) < sum(vclass==2) ;
  vclass = 3-vclass ;
end ;

labels = ones(n,1) ;
labels(~flag) = vclass+1 ;

end

