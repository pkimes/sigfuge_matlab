function SplicePCplotPK2(genename,paramstruct) 
% SplicePCplotPK2, plot for PC loadings building off of splice analysis.
%
%
% NEED TO ALLOW FOR COLORING BY SPLICES AS DEFAULT OPTION AS WELL
%
% ALSO NEED TO ALLOW ANY SORT OF COLORING SCHEME, EVEN SPECIFY MORE THAN 2
%   COLORS
%
% ALSO NEED TO EXTEND SO THAT I CAN HANDLE ANY INPUT OF SAMPLE BAMS AND
% STILL PULL THE NECESSARY FILES (E.G. I SHOULDN'T HAVE TO KNOW WHICH
% SAMPLES I CAN'T ACCESS ETC. -- I SHOULD BE ABLE TO ACCESS THEM ALL...)
%
% 
% STILL NEED TO FIX SO IT CAN TAKE MORE THAN 2 SAMPLES AND CREATE VIOLIN
% PLOTS AND ARROWS FOR THESE COLORS. -- I ALSO NEED TO CHANGE THE ARROW
% BEHAVIOR...
% 
%
% 
% 
% Last updated: 01/31/2013 -- starting to incorporate more splices
%                     1. consider Darshan's 'visually weighted graphs'
%                     2. add lower coverage splices (height by represent?)
%                     3. work w/ raw scale -- what's the benefit of this??
%                     4. Marron's suggestion of plotting at different 
%                        percentiles
%                     5. allowing to change path of ACT graph file
% 
% Last updated: 02/03/2013 -- starting to incorporate more user options
%                     1. users can now specify colors instead of groups
% 
% Last updated: 02/04/2013 -- starting to incorporate more user options
%                     1. allow for vFlag to 1..K
%     NEXT NEED TO ALLOW USER TO REMOVE ARROWS (E.G. WHEN THEY BECOME
%     OVERWHELMING)
%
%
%
% Inputs:
%   genename     - string gene name
%
%   paramstruct  - a Matlab structure of input parameters
%
%     fields           values
%
%      mColors          Kx3 matrix specifying group colors.
%                         error will occur if vFlag ~= 1..size(mColors,2)
%                         default is [1 0 0; 0 0 1] and thus error will
%                         occur if mColors unspecified when 
%                         |max(vFlag)|>2 
%      vFlag           vector of flags for class labels if we don't want
%                         2-means clustering/SigFuge method to determine 
%                         the clusters, default = [], use SigFuge
%      VPidx            length <= 2 vector of which labels to make half
%                         violin plots for,
%                         default is [2,3] when vFlag = []
%                         default is [1,2] when vFlag ~= []
%                         []
%      iVPlot           0,1 indicator of whether to include violin plots,
%                         default = 1
%      iVPuse0          0,1 indicator of whether to use 0's when creating
%                         violin plots (only used if iVPlot = 1),
%                         default = 1
%      iLpit            0,1 indicator of whether to compule L-moments using
%                         prob integral transf to Gaussian inverse CDF,
%                         default = 1
%      iLog             0,1 indicator of whether to use log10 scaled
%                         counts in splice plot, default = 0
%      avgtype          type of average to use when plotting splice plot,
%                         if iLog = 1, average is calculated AFTER log 
%                         transformation, note median will produce same
%                         curves as in coverage depth plots
%                           0 : quantile, default (median) 
%                           1 : mean
%                           2 : winsor mean
%                           3 : trimmed mean
%      avgP             if avgtype == 0, 
%                         quantile to use (1-99), default = 50 (median)
%                       if avgtype == 2 | 3,
%                         percentile to drop from top and bottom when 
%                         calculating winsor or trimmed mean 
%                         (only used if avgtype>=2), i.e. cutoff values at 
%                         [avgP,100-avgP] percentiles, default = 5
%      xUL              define the x upper limit for splice plot
%      minN             minimum number of samples in either cluster for a
%                         splice to be plotted, default = 5
%      minP             minimum proportion of samples in either cluster for a
%                         splice to be plotted, default = 5
%      ACTfile          file name for ACT graph data in 'LUSC_datafiles/ACT/'
%                         path, default = 'ACToutput.txt' 
%      savestr          file name for saving plot, default = genename
%      titlestr         title for plot, default = genename
%
% Output:
%        SplicePlot  - SigFuge curve plot with necessary information
%
% Assumes path can find personal function:
%    vec2matSM.m
%    

%    Copyright (c) P. Kimes 2013

 % general parameters for figure  
mColors = [.3 .3 .3;
           1. .3 .3; ...
           .3 .3 1.];
medColors = max(mColors-.3, 0);
shift = 1;
vFlag = [];
iMedian = 1;
iAll = 1;
iLog = 1;
avgtype = 0;
avgP = 50;
minP = 0.05;
minN = 5;
dir = [];
npc = 4;
minloading = 0.1;
wsplice = 0.5;
plotCurvePC = 1;
plotSplicePC = 1;
depthfile = '';
ACTfile = 'ACToutput.txt';
savestr = genename;
titlestr = [genename ' pile-up curves and splicing data'];
dirstring = 'splice dir';


if nargin > 1;   %  then paramstruct is an argument
    
    if isfield(paramstruct, 'mColors') && ~isempty(paramstruct.mColors);
        mColors = paramstruct.mColors;
        medColors = max(mColors-.4, 0);
    end;

    if isfield(paramstruct, 'shift');
        shift = paramstruct.shift;
    end;

    if isfield(paramstruct, 'iAll');
        iAll = paramstruct.iAll;
    end;
  
    if isfield(paramstruct, 'iMedian');
        iMedian = paramstruct.iMedian;
    end ;

    if isfield(paramstruct, 'vFlag');
        vFlag = paramstruct.vFlag;
    end ;

    if isfield(paramstruct, 'ilog');
        iLog = paramstruct.ilog; 
    end ;

    if isfield(paramstruct, 'avgtype');
        avgtype = paramstruct.avgtype;
        if avgtype ~= 0;
            avgP = 5;
        end;
    end;

    if isfield(paramstruct, 'avgP');
        avgP = paramstruct.avgP;
    end;

    if isfield(paramstruct, 'dir');
        dir = paramstruct.dir; 
    end;

    if isfield(paramstruct, 'wsplice');
        wsplice = paramstruct.wsplice; 
    end;

    if isfield(paramstruct, 'npc');
        npc = paramstruct.npc;
        if npc > 6;
            disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
            disp('!!!  npc must be less than maximum of 6 else !!');
            disp('!!!  plots become too cluttered, using       !!');
            disp('!!!  max of npc = 6                          !!');
            disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
            npc = 6;      
        end;
    end;

    if isfield(paramstruct, 'minloading');
        minloading = paramstruct.minloading; 
        if minloading<0 || minloading>=1;
            disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
            disp('!!!  minloading must be >=0 and <1,          !!');
            disp('!!!  setting to default, minloading = 0.05   !!');
            disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
            minloading = 0.05;
        end;
    end;

    if isfield(paramstruct, 'minN');
        minN = paramstruct.minN;
    end;

    if isfield(paramstruct, 'minP');
        minP = paramstruct.minP;
    end;

    if isfield(paramstruct,'plotCurvePC');
        plotCurvePC = paramstruct.plotCurvePC; 
    end;

    if isfield(paramstruct, 'plotSplicePC');
        plotSplicePC = paramstruct.plotSplicePC; 
    end;

    if isfield(paramstruct, 'depthfile');
        depthfile = paramstruct.depthfile;
    end;

    if isfield(paramstruct, 'ACTfile');
        ACTfile = paramstruct.ACTfile; 
    end;

    if isfield(paramstruct, 'savestr');
        savestr = paramstruct.savestr; 
    end;

    if isfield(paramstruct, 'titlestr');
        titlestr = paramstruct.titlestr; 
    end;

    if isfield(paramstruct, 'dirstring');
        dirstring = paramstruct.dirstring; 
    end;

end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% check whether dimensions of vFlag and mColors match 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(vFlag) ;
  if size(mColors,1)~=length(unique(vFlag)) ;
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
    disp('!!  mColors must be Kx3 color matrix         !!') ;
    disp('!!  where K is the number of unique values   !!') ;
    disp('!!  in vFlag label vector                    !!') ;
    disp('!!  exiting function without any output      !!') ;
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
    return ;
  elseif  ~all(unique(vFlag)==(1:size(mColors,1))) ;
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
    disp('!!  vFlag label vector must consist of values    !!') ;
    disp('!!  between 1..K (K is number of rows in mColors !!') ;
    disp('!!  exiting function without any output          !!') ;
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
    return ;
  end ;
else 
%   VPidx = [2 3] ;
end ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %  the directory where I do everything related to this project
mainpath = ['/Users/pkimes/Dropbox/UNC/Statistics/Research/' ...
            'Projects/NextGen_Splice/MyCode_local/'] ;
% depthpath = [mainpath depthfile] ;
% ACTpath = [mainpath ACTfile] ;
depthpath = depthfile ;
ACTpath = ACTfile ;
exonpath = [mainpath 'LUSC_datafiles/exon_boundaries.txt'] ;
% depthpath = depthfile ;
% ACTpath = ACTfile ;
% exonpath = '/datastore/hayeslab/pkimes/Data/SFdata/exon_boundaries.txt' ;


  %  load in pile-up data
depth = textread(depthpath) ; %#ok
depth = depth' ;

[d,n] = size(depth) ;
logdepth = log10(depth+shift) - log10(shift) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % check to make sure that vFlag is same length as 
  % number of samples
if ~isempty(vFlag) && (length(vFlag) ~= n) ; 
  disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
  disp('!!  vFlag must be nx1 vector.              !!') ;
  disp('!!  exiting function without output.       !!') ;
  disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
  return ;
end ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %  load in splice site data 
  %  (MapSplice 'filtered_normal_junctions.txt' output)
fid = fopen(ACTpath) ;
output = textscan(fid,['%s %s %s %d %d %s',repmat('%f',1,n)]) ;
fclose(fid) ;
tempmat = mat2cell(cell2mat(output(7:end)), ...
                  ones(1,length(output{7})),n) ;
output2 = [output(1:6) {tempmat}] ;
for i = [4 5] ;
  output2(i) = {mat2cell(output2{i},ones(length(output2{i}),1))} ;
end ;
output2 = cat(2,output2{1:end}) ;
colHeadings = {'Chr' 'VType' 'EType' 'Start' 'End' 'Gene' 'Vals'} ;
ACTgraphs = cell2struct(output2,colHeadings,2) ;


%load in exon data
ourAnnot = findinRefSeq(genename, exonpath);

  %  load in exon boundaries 
fid = fopen(exonpath) ;
output = textscan(fid,'%s %s %s %s %s', ...
                      'Delimiter','\t', ...
                      'Headerlines',0, ...
                      'bufsize',100000) ;
fclose(fid) ;
output(2) = {mat2cell(strcmp(output{2},'-')+0,ones(length(output{2}),1))} ;
output(4) = {cellfun(@(x) str2num(x), output{4},'UniformOutput',0)} ; %#ok
output(5) = {cellfun(@(x) str2num(x), output{5},'UniformOutput',0)} ; %#ok
colHeadings = {'Gene' 'NegStrand' 'Chr' 'EStart' 'EEnd'} ;
Annot = cell2struct(cat(2,output{1:end}),colHeadings,2) ;

ourAnnot = Annot(strcmp(genename,{Annot.Gene})) ;
nExons = length(ourAnnot.EStart) ;

  % only parts of structure including the gene that we're interested in
ourACT = ACTgraphs(strcmp(genename,{ACTgraphs.Gene})) ;
ourACTsplices = ourACT(strcmp('splice',{ourACT.EType})) ;


spliceStartExon = cellfun(@(x) ...
                  any(x>=ourAnnot.EStart-5 & x<=ourAnnot.EEnd+5), ...
                num2cell(double([ourACTsplices.Start]))) ;
spliceEndExon = cellfun(@(x) ...
                  any(x>=ourAnnot.EStart-5 & x<=ourAnnot.EEnd+5), ...
                num2cell(double([ourACTsplices.End]))) ;

spliceWeird = ~spliceStartExon & ~spliceEndExon ;


spliceSEsame = cellfun(@(x,y) ...
                  any(x>=ourAnnot.EStart-5 & y<=ourAnnot.EEnd+5), ...
                num2cell(double([ourACTsplices.Start])), ...
                num2cell(double([ourACTsplices.End]))) ;  
              
spliceGoodset = ~spliceWeird & ~spliceSEsame ;

ourACTsplices = ourACTsplices(spliceGoodset) ;

dsplice = length(ourACTsplices) ;
valmat = vertcat(ourACTsplices.Vals) ;

% valmat = valmat(spliceGoodset,:) ;  
% dsplice = sum(spliceGoodset) ;

logvalmat = log10(valmat+shift) - log10(shift) ;    

if isempty(vFlag) ;
    %  obtain cluster labels from 2-means clustering on transformed data
  labels = SigFugeLabelsPK(depth) ;
else
  labels = vFlag ;
end ;
K = max(labels) ;


  %  my guess of the genomic coordinate Chris includes in his pile-up
  %  curves -- not sure why some bases are missing
GenCoord = ourAnnot.EStart(1) ;
for p = 1:nExons ;
    GenCoord = [GenCoord (ourAnnot.EStart(p)+1):ourAnnot.EEnd(p)] ; %#ok
end ;

  %  the exon boundaries for the gene
exons = ourAnnot.EEnd - ourAnnot.EStart ;
exons = [0 cumsum(exons)+1] ;
exonn = zeros(1,d) ;
for p = 1:nExons ;
  if mod(p,2) ;
    exonn((exons(p)+1):exons(p+1)) = 1 ;
  end ;
end ;

  %  need to slightly modify for missing bases in depth coverage files
exons2 = exons ;
exons2(2:(end-1)) = exons2(2:(end-1))-1 ;

  %  add 'introns' between genes to help with visualizing splices
IntronLength = ceil(d/100)*10 ;
plotCurves = [] ;
plotCurves2 = [] ;
plotExonn = [] ;
plotGenCoord = [] ;
iExon = [] ;
for p = 1:nExons ;
  iExon = [iExon ; ...
           ones(exons2(p+1)-exons2(p),1) ; ...
           zeros(IntronLength,1)] ; %#ok
  plotCurves = [plotCurves; ...
                logdepth((exons2(p)+1):exons2(p+1),:); ...
                zeros(IntronLength,n)] ; %#ok
  plotCurves2 = [plotCurves2; ...
                depth((exons2(p)+1):exons2(p+1),:); ...
                zeros(IntronLength,n)] ; %#ok
  plotExonn = [plotExonn ...
               exonn((exons(p)+1):exons(p+1)) ...
               2*ones(1,IntronLength)] ; %#ok
  tempGC1 = (GenCoord(exons(p)+1)-floor(IntronLength/2)): ...
                  (GenCoord(exons(p)+1)-1) ;
  tempGC2 = GenCoord((exons(p)+1):exons(p+1)) ;
  tempGC3 = (GenCoord(exons(p+1))+1): ...
                  (GenCoord(exons(p+1))+ceil(IntronLength/2)) ;
  if p > 1 ;
    tempGC1 = max(tempGC1,ourAnnot.EEnd(p-1)+1) ;
  end ;
  if p < nExons ;
    tempGC3 = min(tempGC3,ourAnnot.EStart(p+1)-1) ;
  end ;    
  plotGenCoord = [plotGenCoord ...
                  tempGC1 ...
                  tempGC2 ...
                  tempGC3 ...
                 ] ; %#ok
end ;
iExon = iExon(1:(end-IntronLength)) ;
plotCurves = plotCurves(1:(end-IntronLength),:) ;
plotCurves2 = plotCurves2(1:(end-IntronLength),:) ;
plotExonn = plotExonn(1:(end-IntronLength)) ;
plotGenCoord = plotGenCoord((floor(IntronLength/2)+1): ...
                            (end-ceil(IntronLength/2))) ;

if ourAnnot.NegStrand ;  %  need to flip if gene is on '-' strand
  exons = exons(end) - fliplr(exons) ;
  plotCurves = flipud(plotCurves) ;
  plotCurves2 = flipud(plotCurves2) ;
  plotExonn = fliplr(plotExonn) ;
  plotGenCoord = fliplr(plotGenCoord) ;
end ;

  %  calculate log medians for all curves
    logmedns = zeros(length(plotGenCoord),K) ;
for k = 1:K ;
  logmedns(:,k) = median(plotCurves(:,labels==k),2) ;
end ;

if iLog ;
  plotCurves2 = log10(plotCurves2+shift) - log10(shift) ;
end ;

  %  use switch cases to define anonymous average function
switch avgtype
  case 0  %  median
    avgFn = @(x) quantile(x,avgP/100,2) ;
  case 1  %  mean
    avgFn = @(x) mean(x,2) ;
  case 2  %  winsor mean
    avgFn = @(x) mean(winsorising(x,100-2*avgP),2) ;
  case 3  %  trimmed mean
    avgFn = @(x) trimmean(x,2*avgP,'round',2) ;
end ;

avg = zeros(length(plotGenCoord),K) ;
for k = 1:K ;
  avg(:,k) = avgFn(plotCurves2(:,labels==k)) ;
end ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % calculate normalized value matricies (i.e. SS = 1)
nlogdepth = plotCurves / sqrt(sum(n*std(logdepth,1,2).^2)) ;
nlogvalmat = logvalmat / sqrt(sum(n*std(logvalmat,1,2).^2)) ;

fullmat = [sqrt(1-wsplice)*nlogdepth ; ...
           sqrt(wsplice)*nlogvalmat] ;

dplot = size(plotCurves,1) ;

  %  defining arrows for splicing information to plot  
if ~isempty(dir) ;
  ddir = size(dir,1) ;
  ndir = size(dir,2) ;
  if ddir ~= (d+dsplice) ;
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
    disp('!!  specified direction vector/matrix    !!') ;
    disp('!!  dimensions differ from data!         !!') ;
    disp('!!  dir matrix must be d+dsplice x ndir. !!') ;
    disp(['!!  d+dsplice = ' num2str(d+dsplice) '               !!']) ;
    disp(['!!  ddir = ' num2str(ddir) '                            !!']) ;
    disp('!!  using default of 4 PC directions     !!') ;
    disp('!!  with curve:splices sum sq. ratio 1:1 !!') ;
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
    dir = [] ;
  else
    npc = ndir ;
    tdir = zeros(dplot,ndir) ;
    tdir(logical(iExon),:) = dir(1:d,:) ;
    if ourAnnot.NegStrand ; 
      tdir = flipud(tdir) ;
    end ;
    dir = [tdir; dir((d+1):(d+dsplice),:)] ;
  end ;
end ;            

if ~isempty(dir) ;
  mdirval = dir ;
else
  paramstruct = struct('npc',npc,'viout',[0 1]) ;
  meigval = pcaSM(fullmat,paramstruct) ;
  mdirval = meigval.meigvec ;
end ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vaxh = axisSM(1:size(plotCurves,1)) ;
vaxv = axisSM(plotCurves) ;  
mmin = min(reshape(plotCurves,1,[])) ;
mmax = max(reshape(plotCurves,1,[])) ;
vaxvB1  = mmin - .06*(mmax - mmin) ;  %  L4 top
vaxvB2  = mmin - .09*(mmax - mmin) ;
vaxvB3  = mmin - .10*(mmax - mmin) ;
% vaxvB4  = mmin - .16*(mmax - mmin) ;
% vaxvB5  = mmin - .19*(mmax - mmin) ;
% vaxvC4  = mmin - .20*(mmax - mmin) ; 
vaxv(1) = -1-2*(npc+1) ;

PCzeros = -1-2*(1:npc) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1) ;
clf ;
hold on ;
axis([vaxh(1) vaxh(2) vaxv(1) vaxv(2)]) ;
ax = axis ;

  %  background rectangle for L4, L6 and exon boundaries
Xeps = (ax(2)-ax(1))*.02 ;

intronCol = hex2dec({'5B' 'A4' 'CB'})' / 255 ;
exonCol = hex2dec({'00' '40' '62'})' / 255 ;

rectangle('Position',[0,vaxvB2,length(plotExonn),vaxvB1-vaxvB2], ...
          'FaceColor',intronCol,'EdgeColor',intronCol) ;
for p = 1:nExons ;
  rectangle('Position',[exons(p)+(p-1)*IntronLength,vaxvB2, ...
                        exons(p+1)-exons(p),vaxvB1-vaxvB2], ...
            'FaceColor',exonCol,'EdgeColor',exonCol) ;
end ;
text(-d/150,(vaxvB1+vaxvB2)/2,ourAnnot.Chr, ...
     'HorizontalAlignment','right', ...
     'VerticalAlignment','middle', ...
     'FontSize',10) ;
    

for i = 1:npc ;
  rectangle('Position',[ax(1)+Xeps,-2*(i+1),(ax(2)-ax(1))-2*Xeps,2], ...
            'FaceColor',[.9 .9 .9]-mod(i,2)*[.2 .2 .2], ...
            'EdgeColor',[.9 .9 .9]-mod(i,2)*[.2 .2 .2], ...
            'curvature',0) ;
end ;

   
  %  adding genomic coordinates
for p = 1:nExons ;
  if ourAnnot.NegStrand ;
    pStart = ourAnnot.EEnd(nExons-p+1) ;
    pEnd = ourAnnot.EStart(nExons-p+1) ;
  else
    pStart = ourAnnot.EStart(p) ;
    pEnd = ourAnnot.EEnd(p) ;
  end ;
  text(exons(p)+(p-1)*IntronLength, ...
       vaxvB3, ...
       num2str(pStart), ...
       'rotation',60, ...
       'horizontalalignment','right', ...
       'verticalalignment','top', ...
       'FontSize',9) ;      
  text(exons(p+1)+(p-1)*IntronLength-20, ...
       vaxvB3, ...
       num2str(pEnd), ...
       'rotation',60, ...
       'horizontalalignment','right', ...
       'verticalalignment','top', ...
       'FontSize',9) ;  
end ;


  %  bottom lines for PC loadings plots
line(vec2matSM([0; length(plotCurves)],npc),vec2matSM(PCzeros,2), ...
     'color','k', ...
     'linewidth',2) ;

   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%%%%%% Calculating the splice plots for PC directions %%%%%%
   
  %  computing splice arrow start and end positions 
  %  p1,p2 = position1 (start) --> position2 (end)
if ourAnnot.NegStrand ;
  p1 = arrayfun(@(x) sum(x<plotGenCoord),[ourACTsplices.End]) ;
  p2 = arrayfun(@(x) sum(x<plotGenCoord),[ourACTsplices.Start]) ;
else
  p1 = arrayfun(@(x) sum(x>plotGenCoord),[ourACTsplices.Start]) ;
  p2 = arrayfun(@(x) sum(x>plotGenCoord),[ourACTsplices.End]) ;        
end ;
p1(p1==0) = -IntronLength/2 ;
  %  extend some splices before gene model
p2(p2==length(plotGenCoord)) = length(plotGenCoord) + IntronLength/2 ;
  %  extend some splice beyond gene model
  
  
  
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %  actually plotting all PC or user specified directions
for i = 1:npc ;

  icolor = (1-mod(i,2))*[.5 0 .5] + mod(i,2)*[0 .5 0] ;
  axis(axis) ;
  if plotSplicePC && max(abs(mdirval((dplot+1):end,i)))>0 ;
    PCarrowH = mdirval((dplot+1):end,i) + PCzeros(i) ;
    vp1 = [p1' PCarrowH] ; vp2 = [p2' PCarrowH] ;
    vp1 = vp1((abs(mdirval((dplot+1):end,i))>=minloading),:) ;
    vp2 = vp2((abs(mdirval((dplot+1):end,i))>=minloading),:) ;
    if ~isempty(vp1) ;
      PCarrows = arrow(vp1,vp2,'Length',4,'TipAngle',30,'Width',1, ...
                      'EdgeColor',icolor,'FaceColor',icolor) ; %#ok
    end ;
  end ;

    % plot the curve loadings (assuming, again, [curves; splices] format
  if plotCurvePC && max(abs(mdirval(1:dplot,i)))>0 ;
    PCcurveH = mdirval(1:dplot, i) / max(abs(mdirval(1:dplot,i))) ...
                + PCzeros(i) ;
    plot(PCcurveH,'-','color',icolor) ;
  end ;

end ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Calculating the splice plots for upper plot %%%%%%

SpliceCnts = zeros(size(valmat,1),K) ;
for k = 1:K ;
  SpliceCnts(:,k) = sum(valmat(:,labels==k)>0,2) ;
end ;

cn = zeros(1,K) ;
for k = 1:K ;
  cn(k) = sum(labels==k) ;
end ;
% 
% if minN > min(cn) ;
%   disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
%   disp('!!  Warning from SplicePCplotPK2.m               !!') ;
%   disp('!!  minN > min(cn), this may result in missing   !!') ;
%   disp('!!  some interesting splices in smallest group   !!') ;
%   disp(['!!  setting minN = min(cn) = ' num2str(min(cn)) '                   !!']) ;
%   disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
%   minN = min(cn) ;
% end ;

fInterest = cellfun(@(x) any(x>=min(floor(cn*minP),minN)), ...
                      num2cell(SpliceCnts,2)) ;
plotACT = ourACTsplices(fInterest) ;

plotvalmat = vertcat(plotACT.Vals) ;
plotvalmat = log10(plotvalmat+shift) - log10(shift) ;
    
mArrowH = zeros(size(plotvalmat,1),K) ;
for k = 1:K ;
  mArrowH(:,k) = avgFn(plotvalmat(:,labels==k)) ;
end ;

  %  add pile-up curves
if iAll ;
  for k = 1:K ;
    plot(plotCurves(:,labels==k),'-','color',mColors(k,:)) ;
  end ;
end ;

  %  add median pile-up curves
if iMedian ;
  for k = 1:K ;
    plot(logmedns(:,k),'-','color',medColors(k,:),'linewidth',2) ;
  end ;
end ;

  %  compute horizontal location of violin plots and dots on splice arrows
xViolins = ((p1(fInterest)+p2(fInterest))/2) ;  %  just centering them b/c I used half violins

vp1 = [p1(fInterest)' mArrowH] ;
vp2 = [p2(fInterest)' mArrowH] ;

  % plotting arrows for splice information
for k = 1:K ;
%   if any(mArrowH(:,k) ~= 0) ;
  axis(axis) ;
  if cn(k)>0 && any(mArrowH(:,k) ~= 0) ;
    A = arrow(vp1(mArrowH(:,k)~=0,[1,k+1]), ...
              vp2(mArrowH(:,k)~=0,[1,k+1]), ...
              'Length',4,'TipAngle',30,'Width',1, ...
              'EdgeColor',medColors(k,:), ...
              'FaceColor',medColors(k,:)) ; %#ok
  end ;
end ;

for k = 1:K ;
  plot(xViolins(mArrowH(:,k)~=0), ...
       mArrowH(mArrowH(:,k)~=0,k), ...
       '.','color',medColors(k,:),'MarkerSize',10+(K-k)*3) ;  %%%%
end ;


  % removing X-axis
set(gca,'XTick',[]) ;
  % re-labeling Y-axis
set(gca,'YTick',[sort(PCzeros) 0:floor(vaxv(2))]) ;
temp = get(gca,'YTickLabel') ;
temp = mat2cell(temp,ones(size(temp,1),1),size(temp,2)) ;
if isempty(dir) ;
  dirtype = 'splice PC' ;
else
  dirtype = dirstring ;
end ;

if length(dirstring) == npc ;
  temp(1:npc) = fliplr(dirstring) ;
else 
  for i = 1:npc ;
    temp(npc+1-i) = {[dirtype num2str(i)]} ;
  end ;
end ;

set(gca,'YTickLabel',temp) ;


  % adding X,Y axis labels and title to plot
xlabel('genomic position') ;
ylabel('log_{10} (counts+1)') ;

if ~isempty(dir) ;
  title(titlestr) ;
else
  title([genename ' pile-up curves and splice PC loadings']) ;
end ;

ax = axis ;
for k = 1:K ;
  text(ax(1)+(ax(2)-ax(1))*.95, ...
       ax(3)+(ax(4)-ax(3))*(1.06-k*.03), ...
       ['cluster ' num2str(k) ' size = ' num2str(cn(k))], ...
       'color',mColors(k,:), ...
       'FontSize',10) ;
end ;


% SaveStr = ['SplicePlot_' savestr] ;
orient landscape ;
% print('-dpdf',SaveStr) ;   
print('-dpdf',savestr) ;   



