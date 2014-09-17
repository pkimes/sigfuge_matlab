function spliceStruct = readinSplices(splicepath, n)
%
%script to pull in splice structure with the following fields:
%   chr,            (string)
%   dir,            (string)
%   gene,           (string)
%   (splice) start, (int)
%   (splice) stop,  (int)
%   vals            (double vector)
%
%written by: Patrick Kimes
%last updated: 10/13/2013


%pull in gene annotations from local file
if nargin < 2;
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
    disp('!!!!!  must specify both splicepath  !!!!!');
    disp('!!!!!  and sample size to work,      !!!!!');
    disp('!!!!!  exiting without output        !!!!!');
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
    return;
end;


%  load in splice site data 
fid = fopen(splicepath);
output = textscan(fid, ...
                  ['%s %s %*s %d %d %s', repmat('%f',1,n)] );
fclose(fid);

nsplice = size(output{1}, 1);

%transform output into structure with each entry corresponding to a splice
tempmat = [output{6:end}];
tempmat = mat2cell(tempmat, ...
                   ones(1, nsplice), ...
                   n);
               
tempbds = num2cell([output{3:4}]);
               
output2 = [output{[1 2 5]} tempbds tempmat];

colHeadings = {'chr' 'dir' 'gene' 'start' 'stop' 'vals'};
spliceStruct = cell2struct(output2, colHeadings, 2);






