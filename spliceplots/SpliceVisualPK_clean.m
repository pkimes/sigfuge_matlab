function SpliceVisualPK_clean(genename,paramstruct) 
% SpliceVisualPK, plot for SigFuge curve view of multiple RNA-seq samples
%                 along gene model.
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
%      iLmom            0, 1 indicator of whether to include L-moment
%                         calculations in figure, default = 1
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
%      iExons           0, 1 indicator of whether to show exon boundaries
%                         or not (just exon 'numbers'), recommended to use
%                         exon #s to prevent overlapping boundary values
%                         when gene model is complex, default = 1
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

%the directory where I do everything related to this project
mainpath = ['/Users/pkimes/Dropbox/UNC/Statistics/Research/' ...
            'Projects/NextGen_Splice/MyCode_local/'] ;

%specify default values
CScluster = 0 ;
mColors = [.3 .3 .3 ;
           1 .3 .3 ; ...
           .3 .3 1] ;
medColors = max(mColors-.3,0) ;
iMedian = 1 ;
vFlag = [] ;
iVPlot = 1 ;
VPidx = [1 2] ;
iVPuse0 = 1 ;
iLmom = 1 ;
iLpit = 1 ;
iLog = 1 ;
avgtype = 0 ;
avgP = 50 ;
xUL = 0 ;
minN = 5 ;
minP = .2 ;
iDiffOnly = 0 ;
cDiffOnly = 1 ;
keyExons = [] ;
iExons = 1 ;
ACTfile = 'ACT/LUSC/ACToutput.txt' ;
depthfile = [mainpath 'LUSC_datafiles/coverage/' genename '_coverage.txt'];
exonfile = [mainpath 'LUSC_datafiles/exon_boundaries.txt'];
savestr = genename ;
titlestr = [genename ' pile-up curves and splicing data'];


if nargin > 1 ;   %  then paramstruct is an argument    
  if isfield(paramstruct,'mColors') && ~isempty(paramstruct.mColors) ;    %  then change to input value
    mColors = paramstruct.mColors ;
    medColors = max(mColors-.3,0) ;
  end ;

  if isfield(paramstruct,'iMedian') ;    %  then change to input value
    iMedian = paramstruct.iMedian ;
  end ;

  if isfield(paramstruct,'CScluster') ;    %  then change to input value
    CScluster = paramstruct.CScluster ;
  end ;

  if isfield(paramstruct,'vFlag') ;    %  then change to input value
    vFlag = paramstruct.vFlag ;
  end ;

  if isfield(paramstruct,'VPidx') ;    %  then change to input value
    VPidx = paramstruct.VPidx ;
    if length(VPidx)>2 ;
      disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
      disp('!! Violin plots can only be made for 2 !!') ;
      disp('!! groups, using only first 2 entries  !!') ;
      disp('!! of VPidx                            !!') ;
      disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
      VPidx = VPidx(1:2) ;
    end ;
  end ;

  if isfield(paramstruct,'iVPlot') ;    %  then change to input value
    iVPlot = paramstruct.iVPlot ;
  end ;

  if isfield(paramstruct,'iVPuse0') ;    %  then change to input value
    iVPuse0 = paramstruct.iVPuse0 ; 
  end ;

  if isfield(paramstruct,'iLmom') ;    %  then change to input value
    iLmom = paramstruct.iLmom ; 
  end ;

  if isfield(paramstruct,'iLpit') ;    %  then change to input value
    iLpit = paramstruct.iLpit ; 
  end ;

  if isfield(paramstruct,'ilog') ;    %  then change to input value
    iLog = paramstruct.ilog ; 
  end ;
  
  if isfield(paramstruct,'avgtype') ;    %  then change to input value
    avgtype = paramstruct.avgtype ;
    if avgtype ~= 0 ;
      avgP = 5 ;
    end ;  
  end ;
  
  if isfield(paramstruct,'avgP') ;    %  then change to input value
    avgP = paramstruct.avgP ; 
  end ;

  if isfield(paramstruct,'xUL') ;    %  then change to input value
    xUL = paramstruct.xUL ; 
  end ;

  if isfield(paramstruct,'minN') ;    %  then change to input value
    minN = paramstruct.minN ; 
  end ;

  if isfield(paramstruct,'minP') ;    %  then change to input value
    minP = paramstruct.minP ; 
  end ;

  if isfield(paramstruct,'iExons') ;    %  then change to input value
    iExons = paramstruct.iExons ; 
  end ;

  if isfield(paramstruct,'iDiffOnly') ;    %  then change to input value
    iDiffOnly = paramstruct.iDiffOnly ; 
  end ;

  if isfield(paramstruct,'cDiffOnly') ;    %  then change to input value
    cDiffOnly = paramstruct.cDiffOnly ; 
  end ;

  if isfield(paramstruct,'keyExons') ;    %  then change to input value
    keyExons = paramstruct.keyExons ; 
  end ;

  if isfield(paramstruct,'ACTfile') ;    %  then change to input value
    ACTfile = paramstruct.ACTfile ; 
  end ;

  if isfield(paramstruct,'depthfile') ;    %  then change to input value
    depthfile = paramstruct.depthfile ; 
  end ;

  if isfield(paramstruct,'exonfile') ;    %  then change to input value
    exonfile = paramstruct.exonfile ; 
  end ;

  if isfield(paramstruct,'savestr') ;    %  then change to input value
    savestr = paramstruct.savestr ; 
  end ;

  if isfield(paramstruct,'titlestr') ;    %  then change to input value
    titlestr = paramstruct.titlestr ; 
  end ;
end ;



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
  VPidx = [2 3] ;
end ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

depthpath = depthfile;
ACTpath = ACTfile;
exonpath = exonfile;


  %  load in pile-up data
depth = textread(depthpath) ; %#ok
depth = depth' ;

[d,n] = size(depth) ;
logdepth = log10(depth+1) ;

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

  % only parts of structure including the gene that we're interested in
ourACT = ACTgraphs(strcmp(genename,{ACTgraphs.Gene})) ;
ourACTsplices = ourACT(strcmp('splice',{ourACT.EType})) ;
ourAnnot = Annot(strcmp(genename,{Annot.Gene})) ;
splice = vertcat(ourACTsplices.Vals) ;

nExons = length(ourAnnot.EStart) ;

if isempty(vFlag) ;
  if CScluster ;
      %  obtain cluster from depth and splice information
      %  50/50 split of total sum sq.
    labels = SigFugeLabelsPK2(depth,splice) ;  
  else
      %  obtain cluster labels from 2-means clustering on transformed data
    labels = SigFugeLabelsPK(depth) ;
  end ;
else
  labels = vFlag ;
end ;
K = max(labels) ;

if K==1 ;
  VPidx = [unique(labels) unique(labels)] ;
end ;

  %  my guess of the genomic coordinate Chris includes in his pile-up
  %  curves -- not sure why some bases are missing
GenCoord = ourAnnot.EStart(1) ;
for p = 1:nExons ;
    GenCoord = [GenCoord (ourAnnot.EStart(p)+1):ourAnnot.EEnd(p)] ; %#ok
end ;

  %  the exon boundaries for the gene
exons = ourAnnot.EEnd - ourAnnot.EStart ;
exons = [0 cumsum(exons)+1] ;

if ~isempty(keyExons) ;
  if min(keyExons)<1 || max(keyExons)>nExons ;
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
    disp('!! keyExons specified exceeds gene model!    !!') ;
    disp('!! figure will be displayed using all exons, !!') ;
    disp('!! re-check exon numbers on figure           !!') ;
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
    keyExons = 1:nExons ;
  end ;
  if ourAnnot.NegStrand ;
    keyExons = nExons+1-keyExons ;
  end ;
  ourAnnot.EStart = ourAnnot.EStart(unique([1 keyExons min(nExons,keyExons+1)])) ;
  ourAnnot.EEnd = ourAnnot.EEnd(unique([keyExons max(1,keyExons-1) nExons])) ;
  exons = exons(unique([1 keyExons keyExons+1 nExons+1])) ;
  nModelExons = nExons ; 
%   nExons = length(ourAnnot.EStart) ;
  nExons = length(exons)-1 ;
end ;

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
for p = 1:nExons ;
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
  plotCurves2 = log10(plotCurves2+1) ;
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
vaxh = axisSM(1:size(plotCurves,1)) ;
vaxv = axisSM(plotCurves) ;  
mmin = min(reshape(plotCurves,1,[])) ;
mmax = max(reshape(plotCurves,1,[])) ;
vaxvB1  = mmin - .06*(mmax - mmin) ;  %  L4 top
vaxvB2  = mmin - .09*(mmax - mmin) ;
vaxvB3  = mmin - .10*(mmax - mmin) ;  %  L4 bottom
vaxvC1  = mmin - .11*(mmax - mmin) ;  %  L6 top
vaxvC2  = mmin - .14*(mmax - mmin) ;
vaxvC3  = mmin - .15*(mmax - mmin) ;  %  L6 bottom
vaxvB4  = mmin - .16*(mmax - mmin) ;
vaxvB5  = mmin - .19*(mmax - mmin) ;
vaxvC4  = mmin - .20*(mmax - mmin) ; 

vaxv(1) = -2 - (vaxv(2)) ;
Splice01 = -1 - ceil(vaxv(2)) ;  %  bottom of splice plot

mmax2 = max(max(avg));  %  max height of average splice curves
sscale = mmax/mmax2 ;  %  scaling factor for Y-axis labels and splice plots
if xUL>0 ;
  sscale = mmax/xUL ;
end ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1) ;
clf ;
hold on ;
axis([vaxh(1) vaxh(2) vaxv(1) vaxv(2)]) ;
ax = axis ;

  %  bottom line for curve plots
line([0 length(plotCurves)],[0 0],'color','k','linewidth',2) ;

  %  background rectangle for L4, L6 and exon boundaries
if iLmom ;
  Xeps = (ax(2)-ax(1))*.02 ;
  Yeps = .1 ;
  rectcolor = [0 139 139] ;
  rectangle('Position',[ax(1)+Xeps,-1+Yeps,(ax(2)-ax(1))-2*Xeps,1-2*Yeps], ...
            'FaceColor',rectcolor./255, ...
            'EdgeColor',rectcolor./255, ...
            'curvature',.4) ;
end ;

  %  add pile-up curves
for k = 1:K ;
  plot(plotCurves(:,labels==k),'-','color',mColors(k,:)) ;
end ;

  %  add median pile-up curves
if iMedian ;
  for k = 1:K ;
    plot(logmedns(:,k),'-','color',medColors(k,:),'linewidth',2) ;
  end ;
end ;

  %  adding exon boundaries to plot
intronCol = hex2dec({'5B' 'A4' 'CB'})' / 255 ;
exonCol = hex2dec({'00' '40' '62'})' / 255 ;
if iLmom ;
  for j = 1:length(plotExonn) ;
    col = exonCol ; % #5BA4CB #004062
    if plotExonn(j) == 2 ;
      col = intronCol ;
    end ;
    plot([j j],[vaxvB4 vaxvB5],'color',col,'linewidth',1) ;
  end ;
  text(-d/175,vaxvB4,ourAnnot.Chr, ...
       'HorizontalAlignment','right', ...
       'VerticalAlignment','top', ...
       'FontSize',10) ;
else
  rectangle('Position',[0, -.08*mmax, length(plotExonn), 0.01*mmax], ...
            'FaceColor', intronCol, 'EdgeColor', intronCol) ;
  for p = 1:nExons ;
    rectangle('Position',[exons(p)+(p-1)*IntronLength, -.10*mmax, ...
                          exons(p+1)-exons(p), .05*mmax], ...
              'FaceColor',exonCol,'EdgeColor',exonCol) ;
  end ;
  text(-d/150,(vaxvB1+vaxvB2)/2,ourAnnot.Chr, ...
       'HorizontalAlignment','right', ...
       'VerticalAlignment','middle', ...
       'FontSize',10) ;
end ;

        
   %  adding L4 values to plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %  calculate L-moments at each base
if iLmom ;
  
  if isempty(vFlag) ;
    if iLpit ;  %  if the prob integral transf L-moments should be used
      vsumstat1 = LstatisticSM(plotCurves(:,labels~=1)',1) ;
      vsumstat2 = LstatisticSM(plotCurves(:,labels~=1)',2) ;
      vsumstat4 = zeros(1,length(vsumstat1)) ;
      vsumstat6 = zeros(1,length(vsumstat1)) ;
      for p = 1:length(vsumstat1) ;
        PITdata = normcdf(plotCurves(p,labels~=1),vsumstat1(p), ...
                            vsumstat2(p)*sqrt(pi)) ;
        vsumstat4(p) = LstatisticSM(PITdata',4) ;
        vsumstat6(p) = LstatisticSM(PITdata',6) ;
      end ;
    else  %  just use the usual L-moments
      vsumstat4 = LstatisticSM(plotCurves(:,labels~=1)',4) ;
      vsumstat6 = LstatisticSM(plotCurves(:,labels~=1)',6) ;    
    end ;
  else
    if iLpit ;  %  if the prob integral transf L-moments should be used
      vsumstat1 = LstatisticSM(plotCurves',1) ;
      vsumstat2 = LstatisticSM(plotCurves',2) ;
      vsumstat4 = zeros(1,length(vsumstat1)) ;
      vsumstat6 = zeros(1,length(vsumstat1)) ;
      for p = 1:length(vsumstat1) ;
        PITdata = normcdf(plotCurves(p,:),vsumstat1(p), ...
                            vsumstat2(p)*sqrt(pi)) ;
        vsumstat4(p) = LstatisticSM(PITdata',4) ;
        vsumstat6(p) = LstatisticSM(PITdata',6) ;
      end ;
    else  %  just use the usual L-moments
      vsumstat4 = LstatisticSM(plotCurves',4) ;
      vsumstat6 = LstatisticSM(plotCurves',6) ;    
    end ;
  end ;

  minv4 = min(vsumstat4) ;
  maxv4 = max(vsumstat4) ;
  for j = 1:length(vsumstat4) ;
    v4j = vsumstat4(j) ;
    icolor = [1/minv4 0 0]*min(v4j,0) + ...
             [0 0 1/maxv4]*max(v4j,0) ;
    icolor = icolor.^(1/2) ;
    plot([j j],[vaxvB1 vaxvB2],'color',icolor,'linewidth',1) ;
  end ;
  tL4 = vsumstat4<0 ;
  for j = 1:length(vsumstat4) ;
    if tL4(j) ;
      plot([j j],[vaxvB2 vaxvB3],'color',[.8 0 .2],'linewidth',1) ;
    else
      plot([j j],[vaxvB2 vaxvB3],'color',[.8 .8 .8],'linewidth',1) ;
    end ;
  end ;
  if iLpit ;
    text(-d/175,vaxvB1,'PIT L4', ...
         'HorizontalAlignment','right','FontSize',9) ;
  else
    text(-d/175,vaxvB1,'L4', ...
         'HorizontalAlignment','right','FontSize',9) ;
  end ;

    %  adding L6 values to plot
  minv6 = min(vsumstat6) ;
  maxv6 = max(vsumstat6) ;
  for j = 1:length(vsumstat4) ;
    v6j = vsumstat6(j) ;
    icolor = [1/minv6 0 0]*min(v6j,0) + ...
             [0 0 1/maxv6]*max(v6j,0) ;
    icolor = icolor.^(1/2) ;
    plot([j j],[vaxvC1 vaxvC2],'color',icolor,'linewidth',1) ;
  end ;
  tL6 = vsumstat6<0 ;
  for j = 1:length(vsumstat4) ;
    if tL6(j) ;
      plot([j j],[vaxvC2 vaxvC3],'color',[.8 0 .2],'linewidth',1) ;
    else
      plot([j j],[vaxvC2 vaxvC3],'color',[.8 .8 .8],'linewidth',1) ;
    end ;
  end ;
  if iLpit ;
    text(-d/175,vaxvC1,'PIT L6', ...
         'HorizontalAlignment','right','FontSize',9) ;
  else
    text(-d/175,vaxvC1,'L6', ...
         'HorizontalAlignment','right','FontSize',9) ;
  end ;
  
end ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %  adding genomic coordinates
if iExons ;
  for p = 1:nExons ;
    if ourAnnot.NegStrand ;
      pStart = ourAnnot.EEnd(nExons-p+1) ;
      pEnd = ourAnnot.EStart(nExons-p+1) ;
    else
      pStart = ourAnnot.EStart(p) ;
      pEnd = ourAnnot.EEnd(p) ;
    end ;
    text(exons(p)+(p-1)*IntronLength, ...
         vaxvC4, ...
         num2str(pStart), ...
         'rotation',60, ...
         'horizontalalignment','right', ...
         'verticalalignment','top', ...
         'FontSize',9) ;      
    text(exons(p+1)+(p-1)*IntronLength-20, ...
         vaxvC4, ...
         num2str(pEnd), ...
         'rotation',60, ...
         'horizontalalignment','right', ...
         'verticalalignment','top', ...
         'FontSize',9) ;  
  end ;
else
  if ~isempty(keyExons) && nExons~=nModelExons ;
    indends = unique([keyExons max(keyExons-1,1) nModelExons]) ;
    txtvals = cell(length(indends),1) ;
    start = 1 ;
    for np = 1:length(indends) ;
      if start == indends(np) ;
        txtvals{np} = ['exon ' num2str(start)] ;
      else
        txtvals{np} = ['exons ' num2str(start) '-' num2str(indends(np))] ;
      end ;
      start = indends(np)+1 ;
    end ;
    for p = 1:nExons ;
      text((exons(p)+exons(p+1))/2+(p-1)*IntronLength, ...
           vaxvC1, ...
           txtvals{p}, ...
           'rotation',90, ...
           'horizontalalignment','right', ...
           'verticalalignment','middle', ...
           'FontSize',9) ;          
    end ;    
  else 
    for p = 1:nExons ;
      text((exons(p)+exons(p+1))/2+(p-1)*IntronLength, ...
           vaxvC1, ...
           ['exon ' num2str(p)], ...
           'rotation',90, ...
           'horizontalalignment','right', ...
           'verticalalignment','middle', ...
           'FontSize',9) ;          
    end ;
  end ;
end ;


  %  bottom line for gene splicing plots
line([0 length(plotCurves)],[Splice01 Splice01],'color','k','linewidth',2) ;

  %  defining arrows for splicing information to plot
% fActualSplice = strcmp('splice',ourAnnot.SEType) ;
% ACTMat = cell2mat(gene.Mat) ;

valmat = vertcat(ourACTsplices.Vals) ;

SpliceCnts = zeros(size(valmat,1),K) ;
for k = 1:K ;
  SpliceCnts(:,k) = sum(valmat(:,labels==k)>0,2) ;
end ;

cn = zeros(1,K) ;
for k = 1:K ;
  cn(k) = sum(labels==k) ;
end ;

fInterest = cellfun(@(x) any(x>=min(floor(cn*minP),minN)), ...
                      num2cell(SpliceCnts,2)) ;
plotACT = ourACTsplices(fInterest) ;


%%%%%% Calculating the splice plots (and violins) for fInterest1 %%%%%%

  %  computing splice arrow start and end positions 
  %  p1,p2 = position1 (start) --> position2 (end)
if ourAnnot.NegStrand ;
  p1 = arrayfun(@(x) sum(x<plotGenCoord),[plotACT.End]) ;
  p2 = arrayfun(@(x) sum(x<plotGenCoord),[plotACT.Start]) ;
else
  p1 = arrayfun(@(x) sum(x>plotGenCoord),[plotACT.Start]) ;
  p2 = arrayfun(@(x) sum(x>plotGenCoord),[plotACT.End]) ;        
end ;
p1(p1==0) = -IntronLength/2 ;
    % extend some splices before gene model
p2(p2==length(plotGenCoord)) = length(plotGenCoord) + IntronLength/2 ;
    % extend some splice beyond gene model
    
plotvalmat = vertcat(plotACT.Vals) ;
if iLog ;  %  if splicing should be plotted on the log10 scale
  plotvalmat = log10(plotvalmat+1) ;
else
  plotvalmat = plotvalmat * sscale ;
end ;

%  shift splice plots down
plotvalmat = plotvalmat + Splice01 ;
    
mArrowH = zeros(size(plotvalmat,1),K) ;
for k = 1:K ;
  mArrowH(:,k) = avgFn(plotvalmat(:,labels==k)) ;
end ;


Violins = cell(sum(fInterest),K) ;
if iVPlot ;  %  whether to include KDE violin plots at all
  for k = VPidx ;
    for j = 1:sum(fInterest) ;
      paramstruct = struct('vh',0, ...
                           'vxgrid',[Splice01; -1], ...
                           'ibdryadj',1, ...
                           'iscreenwrite',0) ;
      tempdata = plotvalmat(j,labels==k) ; 
      if avgtype >= 2 ; % if winsor mean or trimmed mean are used,
                      % drop the bottom bit for violin plots
        cutoff = round(cn(k)*avgP/100) ;
        tempdata = sort(tempdata) ;
        tempdata = tempdata((cutoff+1):(end-cutoff)) ;
      end ;
      if ~iVPuse0 ;  %  whether to include 0s in the KDE violin plots
        tempdata = tempdata(tempdata>Splice01) ;
      end ;
      if length(unique(tempdata)) > 1 ; % may no longer be >1 b/c of no 0s
          [Violins{j,k}(:,2),Violins{j,k}(:,1)] = kdeSM(tempdata',paramstruct) ;
          Violins{j,k}(:,2) = Violins{j,k}(:,2) * length(tempdata) / cn(k) ;
      end ;
    end ;

    Violins(mArrowH(:,k)==Splice01,k) = repmat({[]}, ...
                                          sum(mArrowH(:,k)==Splice01),1) ; % only use violins for arrows
  end ;
end ;

%  compute horizontal location of violin plots and dots on splice arrows
xViolins = ((p1+p2)/2) ;  %  just centering them b/c I used half violins

if ourAnnot.NegStrand ;  %  flip order to avoid error from distributionPlot fn
  Violins = flipud(Violins);
  xViolins = fliplr(xViolins);
  mArrowH = flipud(mArrowH);
  p1 = fliplr(p1);
  p2 = fliplr(p2);
end ;

  % don't show violins for splices that have 'average' coverage really
  % similar between samples
if iDiffOnly ;
  simViolins = abs(mArrowH(:,VPidx(1))-mArrowH(:,VPidx(2))) < cDiffOnly ;
  Violins(simViolins,VPidx) = repmat({[]},sum(simViolins),2) ;
end ;

  % remove violins that correspond to splices not represented in
  % concatenated version of figure
if ~isempty(keyExons) && nExons~=nModelExons ;
  lostsplice = abs(p1-p2) <= 500 ;
  Violins(lostsplice,VPidx) = repmat({[]},sum(lostsplice),2) ;
  mArrowH(lostsplice,:) = Splice01 ;
  disp(['!! number of lost splice = ' num2str(sum(lostsplice))]) ;
  disp(['!! number of showing splices = ' ...
        num2str(sum(mArrowH>Splice01))]) ;
end ;

if all(cn(VPidx)>1) && iVPlot ;
  plotViolins = sum(cellfun(@isempty,Violins(:,VPidx)),2)<2 ;
  Violins = Violins(plotViolins,:) ;
  xViolins = xViolins(plotViolins) ;
  mArrowH = mArrowH(plotViolins,:) ;
  p1 = p1(plotViolins) ;
  p2 = p2(plotViolins) ;
end ;

if iVPlot ;  %  whether to actually include the violin plots
  ax = axis ;
    side = 1 ; sside = 'left' ;
    [xViolinS,sidx] = sort(xViolins) ;
    for k = VPidx ;
      Violins(:,k) = Violins(sidx,k) ;
      eps = (1e-3)*(1:length(Violins(:,k))) ; % in case some violins overlap 
      if ~isempty(Violins(:,k)) ;
        distributionPlot(Violins(:,k), ...
                        'distWidth',.7, ...
                        'variableWidth','false', ...
                        'color',min(mColors(k,:)*3/2+1/4,1), ...
                        'globalNorm',0, ...
                        'groups',[], ...
                        'histOri',sside,...
                        'widthDiv',[2 side],...
                        'histOpt',0, ...
                        'addSpread',0, ...
                        'showMM',0, ...
                        'xValues',xViolinS+eps) ;
      end ;
    side = 2 ; sside = 'right' ;
    end ;
  axis(ax) ;
end ;

checkclass = [2,3] ;
if K < 3 ;
  checkclass = [1,2] ;
end ;
if K > 1 ; 
  for ov = 2:size(mArrowH,1) ;
    overlaps = (p1(1:ov-1)<p2(ov) & p2(1:ov-1)>p1(ov))' ;
    if sum(overlaps)>0 ;
      if mArrowH(ov,checkclass(1)) > Splice01+0.02 ;
        temp = reshape(mArrowH(overlaps,:),[],1);
        temp2 = mArrowH(ov,checkclass(1)) - temp ;
        [delta,idelta] = min(abs(temp2)) ;
        step = 1 ;
        while delta<0.079 && step<5 ;
          if delta == 0 ;
            temp2(idelta) = rand(1)-.5 ;
          end ;
          mArrowH(ov,checkclass(1)) = temp(idelta) + sign(temp2(idelta))*.08 ;
          temp2 = mArrowH(ov,checkclass(1)) - temp ;
          [delta,idelta] = min(abs(temp2)) ;
          step = step+1 ;
        end ;
      end ;
      if mArrowH(ov,checkclass(2)) > Splice01+0.02 ;
        temp = reshape(mArrowH(overlaps,:),[],1);
        temp2 = mArrowH(ov,checkclass(2)) - temp ;
        [delta,idelta] = min(abs(temp2)) ;
        step = 1 ;
        while delta<0.079 && step<5 ;
          if delta == 0 ;
            temp2(idelta) = rand(1)-.5 ;
          end ;
          mArrowH(ov,checkclass(2)) = temp(idelta) + sign(temp2(idelta))*.08 ;
          temp2 = mArrowH(ov,checkclass(2)) - temp ;
          [delta,idelta] = min(abs(temp2)) ;
          step = step+1 ;
        end ;
      end ;    
    end ;
  end ;
end ;

vp1 = [p1' mArrowH] ;
vp2 = [p2' mArrowH] ;

  % plotting arrows for splice information
for k = 2:K ;
  if any(mArrowH(:,k) ~= Splice01) ;
    A = arrow(vp1(mArrowH(:,k)~=(Splice01),[1,k+1]), ...
              vp2(mArrowH(:,k)~=(Splice01),[1,k+1]), ...
              'Length',3,'TipAngle',30,'Width',1, ...
              'EdgeColor',medColors(k,:), ...
              'FaceColor',medColors(k,:)) ; %#ok
  end ;
end ;



  % adding points to connect violin plots with splices
% if ourAnnot.NegStrand ;  %  flip order to avoid error from distributionPlot fn
%   xViolins = fliplr(xViolins) ;
% end ;

for k = 1:K ;
  plot(xViolins(mArrowH(:,k)~=(Splice01)), ...
       mArrowH(mArrowH(:,k)~=(Splice01),k), ...
       '.','color',medColors(k,:),'MarkerSize',10+(K-k)*3) ;  %%%%
end ;

  %  add average curves to the splicing information portion of the plot
if iLog ;
  for k = 1:K ;
    plot(avg(:,k)+Splice01,'-','color',max(0,mColors(k,:)-.1),'linewidth',1);
  end ;  
else 
  for k = 1:K ;
    plot(avg(:,k)*sscale+Splice01,'-','color',max(0,mColors(k,:)-.1),'linewidth',1);
  end ;  
end ;





  %  re-labeling Y-axis and removing X-axis
set(gca,'XTick',[]) ;
yb = Splice01 ;
if iLog ;
  set(gca,'YTick',yb:floor(vaxv(2))) ;
  temp = repmat(num2cell(num2str(0:ceil(vaxv(2)),'%d')),1,2) ;
  set(gca,'YTickLabel',temp) ;
else
  set(gca,'YTick',yb:floor(vaxv(2))) ;
  temp = num2cell([0:1/sscale:(ceil(vaxv(2))/sscale), ...
                  0:ceil(vaxv(2))]) ;
  temp = cellfun(@(x) num2str(x,'%5.0f'),temp,'uniformoutput',0) ;
  set(gca,'YTickLabel',temp) ;
end ;

  %  adding text to plots: X,Y axis labels and Title
xlabel('genomic position') ;
if iLog ;
  ylabel('log_{10} (counts+1)') ;
else
  ylabel('raw counts    /    log_{10} (counts+1)') ;
end ;
title(titlestr,'Interpreter','none') ;

for k = 1:K ;
  text(ax(1)+(ax(2)-ax(1))*.9, ...
       ax(3)+(ax(4)-ax(3))*(1.06-k*.03), ...
       ['cluster ' num2str(k) ' size = ' num2str(cn(k))], ...
       'color',mColors(k,:), ...
       'FontSize',10) ;
end ;
  
SaveStr = ['SplicePlot_' savestr] ;
orient landscape ;
print('-dpdf',SaveStr) ;      



