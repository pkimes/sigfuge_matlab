function SplicePCplotPK(genename,paramstruct) 
% SplicePCplotPK, plot for PC loadings building off of splice analysis.
%
%
% Last updated: 01/28/2013 -- modifying code to allow user to specify
%                             direction vector
% 
% Inputs:
%   genename     - string gene name
%
%   paramstruct  - a Matlab structure of input parameters
%
%     fields           values
%
%      dir              direction matrix of size d x ndir
%      npc              number of principal components to calculate,
%                       default = 4
%      minloading       minimum PC loading (absolute value) to plot arrow,
%                       default = 0.05 
%      iviolins         indicator for to use splice violins and median curves
%                       curves instead of standard bundles of curves,
%                       default = 1
%      iAll             indicator whether to use all splices or just the
%                       subset of interesting ones for calcluating PC
%                       directions, default = 1
%      titlestr         plot title string, default = 
%                         'genename pile-up curves and splice PCs'
%      ACTfile          file name for ACT graph data in 'LUSC_datafiles/ACT/'
%                         path, default = 'ACToutput.txt' 
%      savestr          file name for saving plot, default = genename
%
% Output:
%        SplicePCplot  - SigFuge curve plot with PC loading arrows drawn
%
%
% Assumes path can find personal function:
%    vec2matSM.m
%    

%    Copyright (c) Patrick Kimes 2013

  %  specify default values
dir = [] ;
npc = 4 ;
iviolins = 1 ;
minloading = 0.05 ;
iAll = 1 ;
savestr = genename ;
titlestr = [genename ' pile-up curves and splice loadings'] ;
ACTfile = 'ACToutput.txt' ;

minN = 5 ;
minP = .1 ;

if nargin > 1 ;   %  then paramstruct is an argument      
  if isfield(paramstruct,'npc') ;    %  then change to input value
    npc = paramstruct.npc ; 
    if npc > 6 ;
      disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
      disp('!!!  npc must be less than maximum of 6 else !!') ;
      disp('!!!  plots become too cluttered, using       !!') ;
      disp('!!!  default of npc = 4                      !!') ;
      disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
      npc = 4 ;      
    end ;
  end ;

  if isfield(paramstruct,'iviolins') ;    %  then change to input value
    iviolins = paramstruct.iviolins ; 
    if iviolins ~= 0 && iviolins ~= 1 ;
      disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
      disp('!!!  iviolins indicator must be 0 or 1,      !!') ;
      disp('!!!  setting to default, iviolins = 0        !!') ;
      disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
      iviolins = 0 ;      
    end ;
  end ;

  if isfield(paramstruct,'minloading') ;    %  then change to input value
    minloading = paramstruct.minloading ; 
    if minloading<0 || minloading>=1 ;
      disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
      disp('!!!  minloading must be >=0 and <1,          !!') ;
      disp('!!!  setting to default, minloading = 0.05   !!') ;
      disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
      minloading = 0.05 ;
    end ;
  end ;

  if isfield(paramstruct,'savestr') ;    %  then change to input value
    savestr = paramstruct.savestr ; 
  end ;

  if isfield(paramstruct,'iAll') ;    %  then change to input value
    iAll = paramstruct.iAll ; 
  end ;

  if isfield(paramstruct,'dir') ;    %  then change to input value
    dir = paramstruct.dir ; 
    [ddir,ndir] = size(dir) ;
  end ;

  if isfield(paramstruct,'titlestr') ;    %  then change to input value
    titlestr = paramstruct.titlestr ; 
  end ;

  if isfield(paramstruct,'ACTfile') ;    %  then change to input value
    ACTfile = paramstruct.ACTfile ; 
  end ;

end ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  %  the directory where I do everything related to this project
mainpath = ['/Users/pkimes/Dropbox/UNC/Statistics/Research/' ...
            'Projects/NextGen_Splice/MyCode_local/'] ;

        
  %  load in sample list for splice site data
fid = fopen([mainpath 'LUSC_datafiles/splicesamples178.txt']) ;
SpliceSamples = textscan(fid,'%*s %s','delimiter',',') ;
fclose(fid) ;

SpliceSamples = SpliceSamples{1} ;
fSplice = ones(1,178) ;
fSplice([85 116]) = 0 ; % samples missing according to Darshan
fSplice = logical(fSplice) ;
SpliceSamples = SpliceSamples(fSplice) ;


  %  load in sample list for pile-up data
fid = fopen([mainpath 'LUSC_datafiles/samplenames177.txt']) ;
PileupSamples = textscan(fid,'%s') ;
fclose(fid) ;

PileupSamples = PileupSamples{1} ;
fPileup = ones(1,177) ;
fPileup(115) = 0 ; % samples missing from splice data but in pileup data
fPileup = logical(fPileup) ;
PileupSamples = PileupSamples(fPileup) ;

if sum(strcmp(SpliceSamples,PileupSamples)) == 176 ;
  disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
  disp('!!  splice and pile-up samples match  !!') ;
  disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
else 
  disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
  disp('!!  splice and pile-up samples dont   !!') ;
  disp('!!  match!                            !!') ;
  disp('!!       output is incorrect          !!') ;
  disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
end ;

  %  load in splice site data from Darshan (ACT graph data)
fid = fopen([mainpath 'LUSC_datafiles/ACT/' ACTfile]) ;
output = textscan(fid,['%s %s %s %d %d %s',repmat('%f',1,176)]) ;
fclose(fid) ;

Splice.Chr = output{1} ;
Splice.VType = output{2} ;
Splice.EType = output{3} ;
Splice.Start = output{4} ;
Splice.End = output{5} ;
Splice.Gene = output{6} ;
Splice.Mat = cell2mat(output(7:end)) ;


  %  load in exon boundaries 
fid = fopen([mainpath 'LUSC_datafiles/exon_boundaries.txt']) ;
output = textscan(fid,'%s %s %*s %s %s', ...
                      'Delimiter','\t', ...
                      'Headerlines',0, ...
                      'bufsize',100000) ;
fclose(fid) ;
RefSeq.Gene = output{1} ;
RefSeq.Dir = output{2} ;
RefSeq.Start = output{3} ;
RefSeq.End = output{4} ;
RefSeq.Dir = strcmp(RefSeq.Dir,'-') + 0 ;


  %  starting to look at gene specifically
gene.fSplices = strcmp(genename,Splice.Gene) ;
  gene.SStart = Splice.Start(gene.fSplices) ;
  gene.SEnd = Splice.End(gene.fSplices) ;
  gene.SVType = Splice.VType(gene.fSplices) ;
  gene.SEType = Splice.EType(gene.fSplices) ;
  gene.Mat = mat2cell(Splice.Mat(gene.fSplices,:), ...
                          ones(1,sum(gene.fSplices))) ;
  gene.Chr = unique(Splice.Chr(gene.fSplices)) ;

gene.iRefSeq = find(strcmp(genename,RefSeq.Gene)) ;
  gene.Dir = RefSeq.Dir(gene.iRefSeq) ;
  gene.EStart = str2num(RefSeq.Start{gene.iRefSeq}) ; %#ok<ST2NM>
  gene.EEnd = str2num(RefSeq.End{gene.iRefSeq}) ; %#ok<ST2NM>

gene.nExons = length(gene.EStart) ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  %  load in pile-up data
gene.Curves = textread([mainpath 'LUSC_datafiles/coverage/' ... 
                        genename '_coverage.txt']) ; %#ok
gene.Curves = gene.Curves' ;
gene.Curves = gene.Curves(:,fPileup) ;

[d,n] = size(gene.Curves) ;
gene.logCurves = log10(gene.Curves+1) ;

  %  obtain cluster labels from 2-means clustering on transformed data
medcov = zeros(n,1) ;
percov = zeros(n,1) ;
for i = 1:n ;
  temp = gene.Curves(:,i) ;
  tempflag = (temp == 0) ;
  percov(i) = 1 - sum(tempflag)/d ;
  medcov(i) = median(temp(~tempflag)) ;
end ;
medcov(isnan(medcov)) = 0 ;
flag = (medcov < 5) | (percov < .10) ;
coverage = sum(gene.Curves,1) ;
datanorm = gene.Curves(:,~flag) ./ ...
            vec2matSM(coverage(~flag)+1,d) * median(coverage(~flag)) ;
logdata = log10(datanorm+1) ;

paramstruct = struct('nrep',100, ...
                     'iscreenwrite',0) ;
vclass = SigClust2meanRepSM(logdata,paramstruct) ; 
labSF = logical(vclass-1) ;
if sum(labSF) < sum(~labSF) ;
  labSF = ~labSF ;
end ;

screen = find(~flag) ;
c1 = find(flag) ;
c2 = screen(labSF) ;
c3 = screen(~labSF) ;
% screen = 1:176 ;
% labSF = zeros(1,176) ; labSF(174) = 1 ;
% labSF = logical(labSF) ;
% c1 = [] ;
% c2 = screen(labSF) ;
% c3 = screen(~labSF) ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  %  defining arrows for splicing information to plot
fActualSplice = strcmp('splice',gene.SEType) ;
ACTMat = cell2mat(gene.Mat) ;
valmat = ACTMat(fActualSplice,:) ;
logvalmat = log10(valmat+1) ;

SpliceCount2 = sum(logvalmat(:,c2)>0,2) ;
SpliceCount3 = sum(logvalmat(:,c3)>0,2) ;
fInterest1 = (SpliceCount2 >= min(max(length(c2)*minP,minN),length(c2)) | ...
              SpliceCount3 >= min(max(length(c3)*minP,minN),length(c3))) ;

if iAll ;
  dsplice = sum(fActualSplice) ;
else
  dsplice = sum(fInterest1) ;
end ;
if ~isempty(dir) ;
  if ddir ~= dsplice ;
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
    disp('!!  specified direction vector/matrix  !!') ;
    disp('!!  dimensions differ from data!       !!') ;
    disp('!!  dir matrix must be dsplice x ndir. !!') ;
    disp(['!!  dsplice = ' num2str(dsplice) '                       !!']) ;
    disp(['!!  ddir = ' num2str(ddir) '                          !!']) ;
    disp('!!  using default of 4 PC directions   !!') ;
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
    dir = [] ;
  else
    npc = ndir ;
  end ;
end ;            
if ~isempty(dir) ;
  mdirval = dir ;
  
else
  paramstruct = struct('npc',npc,'viout',[0 1]) ;
  if iAll ;
    meigval = pcaSM(logvalmat,paramstruct) ;
  else
    meigval = pcaSM(logvalmat(fInterest1,:),paramstruct) ;
  end ;
  mdirval = meigval.meigvec ;
end ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  %  my guess of the genomic coordinate Chris includes in his pile-up
  %  curves -- not sure why some bases are missing
GenCoord = gene.EStart(1) ;
for p = 1:gene.nExons ;
    GenCoord = [GenCoord (gene.EStart(p)+1):gene.EEnd(p)] ; %#ok
end ;

  %  the exon boundaries for the gene
exons = gene.EEnd - gene.EStart ;
exons = [0 cumsum(exons)+1] ;
exonn = zeros(1,d) ;
for p = 1:gene.nExons ;
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
plotExonn = [] ;
plotGenCoord = [] ;
for p = 1:gene.nExons ;
  plotCurves = [plotCurves; ...
                gene.logCurves((exons2(p)+1):exons2(p+1),:); ...
                zeros(IntronLength,n)] ; %#ok
  plotExonn = [plotExonn ...
               exonn((exons(p)+1):exons(p+1)) ...
               2*ones(1,IntronLength)] ; %#ok
  tempGC1 = (GenCoord(exons(p)+1)-floor(IntronLength/2)): ...
                  (GenCoord(exons(p)+1)-1) ;
  tempGC2 = GenCoord((exons(p)+1):exons(p+1)) ;
  tempGC3 = (GenCoord(exons(p+1))+1): ...
                  (GenCoord(exons(p+1))+ceil(IntronLength/2)) ;
  if p>1 ;
    tempGC1 = max(tempGC1,gene.EEnd(p-1)+1) ;
  end ;
  if p<gene.nExons ;
    tempGC3 = min(tempGC3,gene.EStart(p+1)-1) ;
  end ;    
  plotGenCoord = [plotGenCoord ...
                  tempGC1 ...
                  tempGC2 ...
                  tempGC3 ...
                 ] ; %#ok
end ;
plotCurves = plotCurves(1:(end-IntronLength),:) ;
plotExonn = plotExonn(1:(end-IntronLength)) ;
plotGenCoord = plotGenCoord((floor(IntronLength/2)+1): ...
                            (end-ceil(IntronLength/2))) ;

if gene.Dir ;  %  need to flip certain things for plotting if gene is on '-' strand
  exons = exons(end) - fliplr(exons) ;
  plotCurves = flipud(plotCurves) ;
  plotExonn = fliplr(plotExonn) ;
  plotGenCoord = fliplr(plotGenCoord) ;
end ;

logmedn1 = median(plotCurves(:,c1),2) ;
logmedn2 = median(plotCurves(:,c2),2) ;
logmedn3 = median(plotCurves(:,c3),2) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vaxh = axisSM(1:size(plotCurves,1)) ;
vaxv = axisSM(plotCurves) ;  
mmin = min(reshape(plotCurves,1,[])) ;
mmax = max(reshape(plotCurves,1,[])) ;
vaxvB4  = mmin - .16*(mmax - mmin) ;
vaxvB5  = mmin - .19*(mmax - mmin) ;
vaxvC4  = mmin - .20*(mmax - mmin) ; 

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
Yeps = .1 ;
rectcolor = [0 139 139] ;
rectangle('Position',[ax(1)+Xeps,-1+Yeps,(ax(2)-ax(1))-2*Xeps,1-2*Yeps], ...
            'FaceColor',rectcolor./255, ...
            'EdgeColor',rectcolor./255, ...
            'curvature',.4) ;

for i = 1:npc ;
  rectangle('Position',[ax(1)+Xeps,-2*(i+1),(ax(2)-ax(1))-2*Xeps,2], ...
            'FaceColor',[.9 .9 .9]-mod(i,2)*[.2 .2 .2], ...
            'EdgeColor',[.9 .9 .9]-mod(i,2)*[.2 .2 .2], ...
            'curvature',0) ;
end ;
          

%  adding exon boundaries to plot
for j = 1:length(plotExonn) ;
  col = [0 .25 .25] ;
  if plotExonn(j) == 2 ;
    col = [127 255 0]./255 ;
  end ;
  plot([j j],[vaxvB4 vaxvB5],'color',col,'linewidth',1) ;
end ;
text(-d/200,vaxvB4,gene.Chr, ...
     'HorizontalAlignment','right', ...
     'VerticalAlignment','top', ...
     'FontSize',10) ;

  %  adding genomic coordinates   
for p = 1:gene.nExons ;
  if gene.Dir ;
    pStart = gene.EEnd(gene.nExons-p+1) ;
    pEnd = gene.EStart(gene.nExons-p+1) ;
  else
    pStart = gene.EStart(p) ;
    pEnd = gene.EEnd(p) ;
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

  %  bottom lines for PC loadings plots
line(vec2matSM([0; length(plotCurves)],npc),vec2matSM(PCzeros,2), ...
     'color','k', ...
     'linewidth',2) ;

   
%  computing splice arrow start and end positions 
%  p1,p2 = position1 (start) --> position2 (end)
tempSEnd = gene.SEnd(fActualSplice) ;
tempSStart = gene.SStart(fActualSplice) ;
if ~iAll ;
  tempSEnd = tempSEnd(fInterest1) ;
  tempSStart = tempSStart(fInterest1) ;
end ;

if gene.Dir ;
  p1 = arrayfun(@(x) sum(x<plotGenCoord),tempSEnd) ;
  p2 = arrayfun(@(x) sum(x<plotGenCoord),tempSStart) ;
else
  p1 = arrayfun(@(x) sum(x>plotGenCoord),tempSStart) ;
  p2 = arrayfun(@(x) sum(x>plotGenCoord),tempSEnd) ;        
end ;
p1(p1==0) = -IntronLength/2 ;
    % extend some splices before gene model
p2(p2==length(plotGenCoord)) = length(plotGenCoord) + IntronLength/2 ;
    % extend some splice beyond gene model

for i = 1:npc ;
    %  shift splice plots down
  PCarrowH = mdirval(:,i) + PCzeros(i) ;

  vp1 = [p1 PCarrowH] ; vp2 = [p2 PCarrowH] ;
  vp1 = vp1((abs(mdirval(:,i))>minloading),:) ;
  vp2 = vp2((abs(mdirval(:,i))>minloading),:) ;

  icolor = (1-mod(i,2))*[.5 0 .5] + mod(i,2)*[0 .5 0] ;
  PCarrows = arrow(vp1,vp2,'Length',4,'TipAngle',30,'Width',1, ...
                  'EdgeColor',icolor,'FaceColor',icolor) ; %#ok
end ;



  %  add main pile-up curves and median curves to the plot
if ~iviolins ;
  plot(plotCurves(:,c1),'-','color',[.4 .4 .4]) ;
  plot(plotCurves(:,c2),'-','color',[1 .3 .3]) ;
  plot(plotCurves(:,c3),'-','color',[.3 .3 1]) ;
  plot(logmedn1,'-','color',[0 0 0],'linewidth',2) ;
  plot(logmedn2,'-','color',[.6 0 0],'linewidth',2) ;
  plot(logmedn3,'-','color',[0 0 .6],'linewidth',2) ;

else 
  
  tempSEnd = gene.SEnd(fActualSplice) ;
  tempSStart = gene.SStart(fActualSplice) ;
  tempMat = gene.Mat(fActualSplice) ;
  
  %  computing splice arrow start and end positions 
  %  p1,p2 = position1 (start) --> position2 (end)
  if gene.Dir ;
    p1 = arrayfun(@(x) sum(x<plotGenCoord),tempSEnd(fInterest1)) ;
    p2 = arrayfun(@(x) sum(x<plotGenCoord),tempSStart(fInterest1)) ;
  else
    p1 = arrayfun(@(x) sum(x>plotGenCoord),tempSStart(fInterest1)) ;
    p2 = arrayfun(@(x) sum(x>plotGenCoord),tempSEnd(fInterest1)) ;        
  end ;
  p1(p1==0) = -IntronLength/2 ;
      % extend some splices before gene model
  p2(p2==length(plotGenCoord)) = length(plotGenCoord) + IntronLength/2 ;
      % extend some splice beyond gene model

  mInterest1 = tempMat(fInterest1) ;
  % mInterest1 is the matrix of splicing values for fInterest1

  c2ArrowH = cellfun(@(x) median(log10(x(c2)+1)),mInterest1) ;
  c3ArrowH = cellfun(@(x) median(log10(x(c3)+1)),mInterest1) ;    

  c2p1 = [p1 c2ArrowH] ; c2p2 = [p2 c2ArrowH] ;
  c3p1 = [p1 c3ArrowH] ; c3p2 = [p2 c3ArrowH] ;


  %  adding violin plots to splicing arrows below curves
  c2vals = cellfun(@(x) log10(x(c2)+1), ...
                  mInterest1,'UniformOutput',0) ;
  c3vals = cellfun(@(x) log10(x(c3)+1), ...
                  mInterest1,'UniformOutput',0) ;

  c2Violins = cell(1,sum(fInterest1)) ;
  c3Violins = cell(1,sum(fInterest1)) ;

  for j = 1:sum(fInterest1) ;
    paramstruct = struct('vh',0, ...
                         'ibdryadj',1, ...
                         'iscreenwrite',0) ;
    datac2 = c2vals{j}' ;
    datac3 = c3vals{j}' ;               
      % note: always >1 if including 0s b/c of "interesting" cutoff
    if length(unique(datac2)) > 1 ; % may no longer be >1 b/c of no 0s
        [c2Violins{j}(:,2),c2Violins{j}(:,1)] = kdeSM(datac2,paramstruct) ;
        c2Violins{j}(:,2) = c2Violins{j}(:,2) * length(datac2) / sum(labSF) ;
    end ;
    if length(unique(datac3)) > 1 ; % may no longer be >1 b/c of no 0s
        [c3Violins{j}(:,2),c3Violins{j}(:,1)] = kdeSM(datac3,paramstruct) ;
        c3Violins{j}(:,2) = c3Violins{j}(:,2) * length(datac3) / sum(~labSF) ;
    end ;
  end ;
  
  c2Violins = c2Violins(c2ArrowH~=0) ; % only use violins for arrows 
  c3Violins = c3Violins(c3ArrowH~=0) ; % above height 0

    %  compute horizontal location of violin plots and dots on splice arrows
  c2xV = ((p1+p2)/2) ;  %  just centering them b/c I used half violins
  c2xV = c2xV(c2ArrowH~=0) ;
  c3xV = ((p1+p2)/2) ;
  c3xV = c3xV(c3ArrowH~=0) ;

  if gene.Dir ;  %  flip order to avoid error from distributionPlot fn
    c2Violins = fliplr(c2Violins) ;
    c3Violins = fliplr(c3Violins) ;
    c2xV = flipud(c2xV) ;  %%%%
    c3xV = flipud(c3xV) ;  %%%%
  end ;

  ax = axis ;
  c2ViolinsT = c2Violins ;
  c3ViolinsT = c3Violins ;
  c2xVT = c2xV ;
  c3xVT = c3xV ;

  [c2xVs,sidx2] = sort(c2xVT) ;
      c2Violins = c2ViolinsT(sidx2) ;
  [c3xVs,sidx3] = sort(c3xVT) ;
      c3Violins = c3ViolinsT(sidx3) ;
  eps2 = (1e-3)*(1:length(c2Violins)) ; % in case some of the violins 
  eps3 = (1e-3)*(1:length(c3Violins)) ; % directly overlap - need space
  distributionPlot(c2Violins, ...
                  'distWidth',.45, ...
                  'variableWidth','false', ...
                  'color',[1 .4 .4], ...
                  'globalNorm',0, ...
                  'groups',[], ...
                  'histOri','left',...
                  'widthDiv',[2 1],...
                  'histOpt',0, ...
                  'addSpread',0, ...
                  'showMM',0, ...
                  'xValues',c2xVs+eps2') ;
  distributionPlot(c3Violins, ...
                  'distWidth',.35, ...
                  'variableWidth','false', ...
                  'color',[.4 .4 1], ...
                  'globalNorm',0, ...
                  'groups',[], ...
                  'histOri','right',...
                  'widthDiv',[2 2],...
                  'histOpt',0, ...
                  'addSpread',0, ...
                  'showMM',0, ...
                  'xValues',c3xVs+eps3') ;
  axis(ax) ;

%   [A,B] = meshgrid(c2ArrowH,c3ArrowH) ;

  % plotting arrows for splice information
  c2p1 = c2p1(c2ArrowH~=0,:) ; c2p2 = c2p2(c2ArrowH~=0,:) ;  %%%%
  c3p1 = c3p1(c3ArrowH~=0,:) ; c3p2 = c3p2(c3ArrowH~=0,:) ;  %%%%

  c2Arrows = arrow(c2p1,c2p2,'Length',4,'TipAngle',30,'Width',1, ...
                  'EdgeColor',[.6 0 0],'FaceColor',[.6 0 0]) ; %#ok
  c3Arrows = arrow(c3p1,c3p2,'Length',4,'TipAngle',30,'Width',1, ...
                  'EdgeColor',[0 0 .6],'FaceColor',[0 0 .6]) ; %#ok

    % adding points to connect violin plots with splices
  if gene.Dir ;  %  flip order to avoid error from distributionPlot fn
    c2xV = flipud(c2xV) ;
    c3xV = flipud(c3xV) ;
  end ;
  % plot(c2xV,c2ArrowH,'.','color',[.6 0 0],'MarkerSize',20) ;
  % plot(c3xV,c3ArrowH,'.','color',[0 0 .6],'MarkerSize',20) ;
  plot(c2xV,c2ArrowH(c2ArrowH~=0),'.','color',[.6 0 0],'MarkerSize',20) ;  %%%%
  plot(c3xV,c3ArrowH(c3ArrowH~=0),'.','color',[0 0 .6],'MarkerSize',15) ;  %%%%

  plot(logmedn1,'-','color',[.4 .4 .4],'linewidth',2) ;
  plot(logmedn2,'-','color',[.8 .4 .4],'linewidth',2) ;
  plot(logmedn3,'-','color',[.4 .4 .8],'linewidth',2) ;

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
  dirtype = 'splice dir' ;
end ;

for i = 1:npc ;
  temp(npc+1-i) = {[dirtype num2str(i)]} ;
end ;
set(gca,'YTickLabel',temp) ;


  %  adding text to plots: min/max L4/L6 and sample sizes, 
  %  also adding X,Y axis labels and Title
xlabel('genomic position') ;
ylabel('log_{10} (counts+1)') ;

if ~isempty(dir) ;
  title(titlestr) ;
else
  title([genename ' pile-up curves and splice PC loadings']) ;
end ;

text(ax(1)+(ax(2)-ax(1))*.9, ...
     ax(3)+(ax(4)-ax(3))*1.06, ...
     '#(Red,Blue,Black)/Total:', ...
     'FontSize',10) ;
text(ax(1)+(ax(2)-ax(1))*.9, ...
     ax(3)+(ax(4)-ax(3))*1.03, ...
     ['(' num2str(length(c2)) ',' num2str(length(c3)) ',' ...
     num2str(length(c1)) ') / ' num2str(n)], ...
      'FontSize',10) ;              


SaveStr = ['GraphicalOutput/SplicePlot_' savestr] ;
orient landscape ;
print('-dpdf',SaveStr) ;      



