function filterset = SpliceFilterPK(splices, annots)
%SPLICEFILTERPK takes in splice structure and exon annotation pair
%               for a single gene - returns filtered set of splices
%
% Inputs:
%   splices  - a Matlab structure containing splicing information,
%              typically output form findinSplices.m
%
%   annots   - a Matlab structure containing exon boundary information,
%              typically output from findinRefSeq.m
%
% Output:
%   filterstruct  - splices structure filtered to smaller subset
%
%
%written by: Patrick Kimes
%last updated: 10/13/2013


%current filtering procedure:

%1. remove splices that don't originate and end
%   within or within 5bp of the exon boundaries
startinexon = cellfun(@(x) ...
                  any(x>=annots.start-5 & x<=annots.stop+5), ...
                num2cell(double([splices.start])));
stopinexon = cellfun(@(x) ...
                  any(x>=annots.start-5 & x<=annots.stop+5), ...
                num2cell(double([splices.stop])));

spliceOut = ~startinexon & ~stopinexon;


%2. remove splices that originate and end in the sampe exon
%
sameStartStop = cellfun(@(x,y) ...
                  any(x>=annots.start-5 & y<=annots.stop+5), ...
                num2cell(double([splices.start])), ...
                num2cell(double([splices.stop])));  
              

%combine the rules
spliceGoodset = ~spliceOut & ~sameStartStop;


%subset the splices
filterset = splices(spliceGoodset);

          
          
end

