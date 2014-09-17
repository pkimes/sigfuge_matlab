function SpliceVisualPK2(genename,paramstruct) 
% SpliceVisualPK2, plot for curve view of multiple RNA-seq samples
%                  along gene model.
%
%
% Inputs:
%   genename     - string gene name
%
%   paramstruct  - a Matlab structure of input parameters
%
%     fields           values
%      iVPlot           0,1 indicator of whether to include violin plots,
%                         default = 1
%      VPuse0           0,1 indicator of whether to use 0's when creating
%                         violin plots (only used if iVPlot = 1),
%                         default = 1
%      iLpit            0,1 indicator of whether to compule L-moments using
%                         prob integral transf to Gaussian inverse CDF,
%                         default = 1
%      ilog             0,1 indicator of whether to use log10 scaled
%                         counts in splice plot, default = 0
%      avgtype          type of average to use when plotting splice plot,
%                         if ilog = 1, average is calculated AFTER log 
%                         transformation, note median will produce same
%                         curves as in coverage depth plots
%                           0 : median, default 
%                           1 : mean
%                           2 : winsor mean
%                           3 : trimmed mean
%      avgP             percentile to drop when calculating winsor or
%                         trimmed mean (only used if avgtype>=2), i.e.
%                         cutoff values at [avgP,100-avgP] percentiles,
%                         default = 5
%      xUL              define the x upper limit for splice plot
%      minN             minimum number of samples in each cluster for a
%                         splice to be plotted, default = 5
%      savestr          file name for saving plot, default = genename
%
% Output:
%        SplicePlot  - SigFuge curve plot with necessary information
%
% Assumes path can find personal function:
%    vec2matSM.m
%    

%    Copyright (c) P. Kimes 2012

  %  specify default values
iVPlot = 1 ;
iVPuse0 = 1 ;
iLpit = 1 ;
ilog = 0 ;
avgtype = 0 ;
avgP = 5 ;
xUL = 0 ;
minN = 5 ;
savestr = genename ;

if nargin > 1 ;   %  then paramstruct is an argument    
  if isfield(paramstruct,'iVPlot') ;    %  then change to input value
    iVPlot = paramstruct.iVPlot ;
  end ;

  if isfield(paramstruct,'iVPuse0') ;    %  then change to input value
    iVPuse0 = paramstruct.iVPuse0 ; 
  end ;

  if isfield(paramstruct,'iLpit') ;    %  then change to input value
    iLpit = paramstruct.iLpit ; 
  end ;

  if isfield(paramstruct,'ilog') ;    %  then change to input value
    ilog = paramstruct.ilog ; 
  end ;
  
  if isfield(paramstruct,'avgtype') ;    %  then change to input value
    avgtype = paramstruct.avgtype ; 
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

  if isfield(paramstruct,'savestr') ;    %  then change to input value
    savestr = paramstruct.savestr ; 
  end ;
end ;


mainpath = ['/Users/pkimes/Dropbox/UNC/Statistics/Research/' ...
            'Projects/NextGen_Splice/MyCode_local/'] ;

        
  %  load in sample list for splice site data
fid = fopen([mainpath 'LUSC_datafiles/splicesamples178.txt']) ;
SpliceSamples = textscan(fid,'%*s %s','delimiter',',') ;
fclose(fid) ;

SpliceSamples = SpliceSamples{1} ;
fSplice = zeros(1,178) ;
fSplice([85 116]) = 1 ; % samples missing according to Darshan
fSplice = logical(fSplice) ;
SpliceSamples = SpliceSamples(~fSplice) ;


  %  load in sample list for pile-up data
fid = fopen([mainpath 'LUSC_datafiles/samplenames177.txt']) ;
PileupSamples = textscan(fid,'%s') ;
fclose(fid) ;

PileupSamples = PileupSamples{1} ;
fPileup = zeros(1,177) ;
fPileup(115) = 1 ; % samples missing from splice data but in pileup data
fPileup = logical(fPileup) ;
PileupSamples = PileupSamples(~fPileup) ;

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
fid = fopen([mainpath 'LUSC_datafiles/ACToutput.txt']) ;
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


  %  load in pile-up data
gene.Curves = textread([mainpath 'LUSC_datafiles/coverage/' ... 
                        genename '_coverage.txt']) ; %#ok
gene.Curves = gene.Curves' ;
gene.Curves = gene.Curves(:,~fPileup) ;

[d,n] = size(gene.Curves) ;
gene.logCurves = log10(gene.Curves+1) ;


  %  obtain cluster labels from 2-means clustering on transformed data
% medcov = zeros(n,1) ;
% percov = zeros(n,1) ;
% for i = 1:n ;
%     temp = gene.Curves(:,i) ;
%     tempflag = (temp == 0) ;
%     percov(i) = 1 - sum(tempflag)/d ;
%     medcov(i) = median(temp(~tempflag)) ;
% end ;
% medcov(isnan(medcov)) = 0 ;
% flag = (medcov < 5) | (percov < .10) ;
% coverage = sum(gene.Curves,1) ;
% datanorm = gene.Curves(:,~flag) ./ ...
%             vec2matSM(coverage(~flag)+1,d) * median(coverage(~flag)) ;
% logdata = log10(datanorm+1) ;
% 
% paramstruct = struct('nrep',100, ...
%                      'iscreenwrite',0) ;
% vclass = SigClust2meanRepSM(logdata,paramstruct) ; 
% labSF = logical(vclass-1) ;
% if sum(labSF)/sum(~flag) < .5 ;
%     labSF = ~labSF ;
% end ;


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
plotCurves2 = [] ;
plotExonn = [] ;
plotGenCoord = [] ;
for p = 1:gene.nExons ;
    plotCurves = [plotCurves; ...
                  gene.logCurves((exons2(p)+1):exons2(p+1),:); ...
                  zeros(IntronLength,n)] ; %#ok
    plotCurves2 = [plotCurves2; ...
                  gene.Curves((exons2(p)+1):exons2(p+1),:); ...
                  zeros(IntronLength,n)] ; %#ok
    plotExonn = [plotExonn ...
                 exonn((exons(p)+1):exons(p+1)) ...
                 2*ones(1,IntronLength)] ; %#ok
    plotGenCoord = [plotGenCoord ...
                    (GenCoord(exons(p)+1)-floor(IntronLength/2)): ...
                        (GenCoord(exons(p)+1)-1) ...
                    GenCoord((exons(p)+1):exons(p+1)) ...
                    (GenCoord(exons(p+1))+1): ...
                        (GenCoord(exons(p+1))+ceil(IntronLength/2)) ...
                   ] ; %#ok
end ;
plotCurves = plotCurves(1:(end-IntronLength),:) ;
plotCurves2 = plotCurves2(1:(end-IntronLength),:) ;
plotExonn = plotExonn(1:(end-IntronLength)) ;
plotGenCoord = plotGenCoord((floor(IntronLength/2)+1): ...
                            (end-ceil(IntronLength/2))) ;

if gene.Dir ;  %  need to flip certain things for plotting if gene is on '-' strand
    exons = exons(end) - fliplr(exons) ;
    plotCurves = flipud(plotCurves) ;
    plotCurves2 = flipud(plotCurves2) ;
    plotExonn = fliplr(plotExonn) ;
    plotGenCoord = fliplr(plotGenCoord) ;
end ;

% screen = find(~flag) ;
% c1 = find(flag) ;
% c2 = screen(labSF) ;
% c3 = screen(~labSF) ;
% logmedn1 = median(plotCurves(:,c1),2) ;
% logmedn2 = median(plotCurves(:,c2),2) ;
% logmedn3 = median(plotCurves(:,c3),2) ;
logmedn = median(plotCurves,2) ;

if ilog ;
    plotCurves2 = log10(plotCurves2+1) ;
end ;
    
if avgtype == 0 ;  %  use median for splice plots
    avg = median(plotCurves2,2) ;
elseif avgtype == 1 ;  %  use mean for splice plots
    avg = mean(plotCurves2,2) ;
elseif avgtype == 2 ;  %  use winsor mean for splice plots
    avg = mean(winsorising(plotCurves2,100-2*avgP),2) ;
elseif avgtype == 3 ;  %  use trimmed mean for splice plots
    avg = trimmean(plotCurves2,avgP,'round',2) ;
end ;

if iLpit;  %  if the prob integral transf L-moments should be used
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vaxh = axisSM(1:size(plotCurves,1)) ;
vaxv = axisSM(plotCurves) ;  
mmin = min(reshape(plotCurves,1,[])) ;
mmax = max(reshape(plotCurves,1,[])) ;
vaxvB1  = mmin - .06*(mmax - mmin) ;  %  L4 top
vaxvB2  = mmin - .09*(mmax - mmin) ;
vaxvB3  = mmin - .10*(mmax - mmin) ;  %  L4 bottom
vaxvB4  = mmin - .16*(mmax - mmin) ;
vaxvB5  = mmin - .19*(mmax - mmin) ;
vaxvC1  = mmin - .11*(mmax - mmin) ;  %  L6 top
vaxvC2  = mmin - .14*(mmax - mmin) ;
vaxvC3  = mmin - .15*(mmax - mmin) ;  %  L6 bottom
vaxvC4  = mmin - .20*(mmax - mmin) ; 
vaxv(1) = -2 - (vaxv(2)) ;
Splice0 = -1 - ceil(vaxv(2)) ;
mmax2 = max(avg);  %  max height of average splice curves
sscale = mmax/mmax2 ;  %  scaling factor if ilog = 0
if xUL>0 ;
    sscale = mmax/xUL ;
end ;

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
        

  %  add main pile-up curves and median curves to the plot
plot(plotCurves,'-') ;
plot(logmedn,'-','color',[.2 .2 .2],'linewidth',3) ;


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


   %  adding L4 values to plot
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
    text(-d/200,vaxvB1,'PIT L4', ...
         'HorizontalAlignment','right','FontSize',9) ;
else
    text(-d/200,vaxvB1,'L4', ...
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
    text(-d/200,vaxvC1,'PIT L6', ...
         'HorizontalAlignment','right','FontSize',9) ;
else
    text(-d/200,vaxvC1,'L6', ...
         'HorizontalAlignment','right','FontSize',9) ;
end ;


for p = 1:gene.nExons ;  %  adding genomic coordinates
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


  %  bottom line for gene splicing plots
line([0 length(plotCurves)],[Splice0 Splice0],'color','k','linewidth',2) ;


  %  defining arrows for splicing information to plot
fActualSplice = strcmp('splice',gene.SEType) ;
SpliceMat = cell2mat(gene.Mat) ;
SpliceCount = sum(SpliceMat>0,2) ;
fInterest = (SpliceCount>=minN) & fActualSplice ;

  %  computing splice arrow start and end positions 
if gene.Dir ;
    p1 = arrayfun(@(x) sum(x<plotGenCoord),gene.SEnd(fInterest)) ;
    p2 = arrayfun(@(x) sum(x<plotGenCoord),gene.SStart(fInterest)) ;
else
    p1 = arrayfun(@(x) sum(x>plotGenCoord),gene.SStart(fInterest)) ;
    p2 = arrayfun(@(x) sum(x>plotGenCoord),gene.SEnd(fInterest)) ;        
end ;
p1(p1==0) = -IntronLength/2 ;
p2(p2==length(plotGenCoord)) = length(plotGenCoord) + IntronLength/2 ;

if ilog ;  %  if splicing should be plotted on the log10 scale
    mInterest = cellfun(@(x) log10(x+1),gene.Mat(fInterest),'UniformOutput',0) ;
else
    mInterest = gene.Mat(fInterest) ;
end ;
%     c2ArrowH = cellfun(@(x) median(log10(x(c2)+1)),gene.Mat(fInterest))+Splice0 ;
%     c3ArrowH = cellfun(@(x) median(log10(x(c3)+1)),gene.Mat(fInterest))+Splice0 ;
% else
if avgtype == 0 ;  %  use median for splice plots
    ArrowH = cellfun(@(x) median(x),mInterest) ;
elseif avgtype == 1 ;  %  use mean for splice plots
    ArrowH = cellfun(@(x) mean(x),mInterest) ;
elseif avgtype == 2 ;  %  use winsor mean for splice plots
    ArrowH = cellfun(@(x) mean(winsorising(x,100-2*avgP)),mInterest) ;
elseif avgtype == 3 ;  %  use trimmed mean for splice plots
    ArrowH = cellfun(@(x) trimmean(x,avgP,'round'),mInterest) ;
end ;
% end ;

if ~ilog ;  %  no need to scale if on log 
    ArrowH = ArrowH * sscale ;
end ;

  %  shift splice plots down
ArrowH = ArrowH + Splice0 ;

vp1 = [p1 ArrowH] ; vp2 = [p2 ArrowH] ;


  %  adding violin plots to splicing arrows below curves
if ilog ;  %  if splicing should be plotted on the log10 or raw scale
    vvals = cellfun(@(x) log10(x+1)+Splice0, ...
                    gene.Mat(fInterest),'UniformOutput',0) ;
else
    vvals = cellfun(@(x) x*sscale+Splice0, ...
                    gene.Mat(fInterest),'UniformOutput',0) ;
end ;

vViolins = cell(1,sum(fInterest)) ;

if iVPlot ;  %  whether to include KDE violin plots at all
    for j = 1:sum(fInterest) ;
        paramstruct = struct('vh',0, ...
                             'ibdryadj',0, ...
                             'iscreenwrite',0) ;
        if ~iVPuse0 ;  %  whether to include 0s in the KDE violin plots
            vdata = vvals{j}' ; vdata = vdata(vdata>Splice0) ;
            [vViolins{j}(:,2),vViolins{j}(:,1)] = kdeSM(vdata,paramstruct) ;
            vViolins{j}(:,2) = vViolins{j}(:,2) * length(vdata) / n ;
        else
            [vViolins{j}(:,2),vViolins{j}(:,1)] = kdeSM(vvals{j}',paramstruct) ;
        end ;    
    end ;

%     c2Violins = c2Violins(c2ArrowH~=(Splice0)) ;
%     c3Violins = c3Violins(c3ArrowH~=(Splice0)) ;
end ;

  %  compute horizontal location of violin plots and dots on splice arrows
vxV = ((p1+p2)/2) ; 

if gene.Dir ;  %  flip order to avoid error from distributionPlot fn
    vViolins = fliplr(vViolins) ;
    vxV = flipud(vxV) ;
end ;

if iVPlot ;  %  whether to actually include the violin plots
    ax = axis ;
    if ~iVPuse0 ;  %  whether to use 0s in KDE volin plots
        mcolor = repmat({[.5 .5 .5]},1,length(vViolins)) ;
        [~,sidx] = sort(vxV) ;
        mcolorT = mcolor(sidx) ;
        vViolinsT = vViolins(sidx) ;
        vxVT = vxV(sidx) ;
        distributionPlot(vViolinsT, ...
                        'distWidth',5, ...
                        'variableWidth','false', ...
                        'color',mcolorT, ...
                        'globalNorm',1, ...
                        'groups',[], ...
                        'histOpt',0, ...
                        'divFactor',[25,2,1], ...
                        'addSpread',0, ...
                        'showMM',0, ...
                        'xValues',vxVT) ;
    else
        [vxVT,sidx] = sort(vxV) ;
        vViolinsT = vViolins(sidx) ;
        eps = (1e-3)*(1:length(vViolins)) ;
        distributionPlot(vViolinsT, ...
                        'distWidth',10, ...
                        'variableWidth','false', ...
                        'color',[.5 .5 .5], ...
                        'globalNorm',0, ...
                        'groups',[], ...
                        'histOpt',0, ...
                        'divFactor',[25,2,1], ...
                        'addSpread',0, ...
                        'showMM',0, ...
                        'xValues',vxVT+eps') ;

    end ;
    axis(ax) ;
end ;

  % plotting arrows for splice information
% c2p1 = c2p1(c2ArrowH~=(Splice0),:) ; c2p2 = c2p2(c2ArrowH~=(Splice0),:) ;
% c3p1 = c3p1(c3ArrowH~=(Splice0),:) ; c3p2 = c3p2(c3ArrowH~=(Splice0),:) ;

vArrows = arrow(vp1,vp2,'Length',4,'TipAngle',30,'Width',1, ...
                'EdgeColor',[.3 .3 .3],'FaceColor',[.3 .3 .3]) ; %#ok

  % adding points to connect violin plots with splices
% plot(flipud(c2xV),c2ArrowH(c2ArrowH~=(Splice0)),'.','color',[.6 0 0],'MarkerSize',20) ;
% plot(flipud(c3xV),c3ArrowH(c3ArrowH~=(Splice0)),'.','color',[0 0 .6],'MarkerSize',20) ;
if gene.Dir ;  %  flip order to avoid error from distributionPlot fn
    vxV = flipud(vxV) ;
end ;
plot(vxV,ArrowH,'.','color',[.3 .3 .3],'MarkerSize',20) ;

  %  add main pile-up curves and median curves to the plot
if ilog ;
    plot(avg+Splice0,'-','color',[.4 .4 .4],'linewidth',1) ;
else 
    plot(avg*sscale+Splice0,'-','color',[.4 .4 .4],'linewidth',1) ;
end ;

  %  re-labeling Y-axis and removing X-axis
set(gca,'XTick',[]) ;
if ilog ;
    set(gca,'YTick',Splice0:floor(vaxv(2))) ;
    set(gca,'YTickLabel', ...
        repmat(num2cell(num2str(0:ceil(vaxv(2)),'%d')),1,2)) ;
else
    set(gca,'YTick',Splice0:floor(vaxv(2))) ;
    temp = num2cell([0:1/sscale:(ceil(vaxv(2))/sscale),0:ceil(vaxv(2))]) ;
    temp = cellfun(@(x) num2str(x,'%5.0f'),temp,'uniformoutput',0) ;
    set(gca,'YTickLabel',temp) ;
end ;

  %  adding text to plots: min/max L4/L6 and sample sizes, 
  %  also adding X,Y axis labels and Title
xlabel('genomic position') ;
if ilog ;
    ylabel('log_{10} (counts+1)') ;
else
    ylabel('raw counts    /    log_{10} (counts+1)') ;
end ;
title([genename ' pile-up curves and splicing data']) ;
% text(ax(1)+(ax(2)-ax(1))*.9, ...
%      ax(3)+(ax(4)-ax(3))*1.06, ...
%      '#(Red,Blue,Black)/Total:', ...
%      'FontSize',10) ;
% text(ax(1)+(ax(2)-ax(1))*.9, ...
%      ax(3)+(ax(4)-ax(3))*1.03, ...
%      ['(' num2str(length(c2)) ',' num2str(length(c3)) ',' ...
%      num2str(length(c1)) ') / ' num2str(n)], ...
%       'FontSize',10) ;              

  
SaveStr = ['GraphicalOutput/SplicePlot2_' savestr] ;
orient landscape ;
print('-dpdf',SaveStr) ;      




