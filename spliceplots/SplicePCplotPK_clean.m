function SplicePCplotPK_clean(genename, ...
                              depthpath, ...
                              splicepath, ...
                              exonpath, ...
                              paramstruct) 
% SplicePCplotPK_clean, plot for PC loadings combining splicing and
% coverage information
%
%
% Inputs:
%   genename     - string gene name
%   splicepath   - file path for splice data, default = '' 
%   depthpath    - file path for coverage data, default = ''
%   exonpath     - file path for exon path, default = ''
%   paramstruct  - a Matlab structure of input parameters
%      -------------    ----------------------
%      fields           values
%      -------------    ----------------------
%  general plotting
%      mColors          Kx3 matrix specifying group colors,
%                         default is equally spaced colors on HSV colormap
%                         except if K=3, default is: black,red,blue
%      labels           vector of class labels if we don't want
%                         2-means clustering/SigFuge method to determine 
%                         the clusters, default = [], will use SigFuge
%      savestr          file name for saving plot, default = genename
%      titlestr         title for plot, default = genename
%
%  calculations
%      iLog             0,1 indicator of whether to use log10 scaled
%                         counts in splice plot, default = 0
%      shift            numeric value to use as offset in log
%                         transformations, i.e. log(x+shift) - log(shift),
%                         default = 1
%      iAvg             0,1 indicator of whether to include average curves
%                         in top panel, default = 1
%      iAll             0,1 indicator of whether to include all n curves
%                         in top panel, default = 0
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
%      iexbds           0,1 indicator whether to include chromosome
%                         position numbers at exon boundaries, default = 1
%
%  bottom dir plots
%      dir              numeric matrix specifying directions to be plotted,
%                         dir should be a ddir x ndir matrix where 
%                         ddir = (length of curves) + (# of splices),
%                         default = [] will calculate npc dirs by PCA
%      npc              numeric value specifying number of PC directions to
%                         calculate, has no effect when dir ~= [], value is
%                         upperbounded at 6, default = 4
%      minloading       numeric value in [0,1) for minimum abs loading of a
%                         splice to be plotted in the direction plots, used
%                         to reduce clutter w/ a lot of small splices,
%                         default = 0.1
%      wsplice          numeric value in [0,1] for % of total sum sq. that
%                         should be allocated to splices, 
%                         e.g. 0 corresponds to PCA on only curve data,
%                         and 1 corresponds to PCA on only splicing data,
%                         default = 0.5 (equal)
%      plotCurvePC      0,1 indicator of whether to show curves in
%                         direction plots, default = 1
%                         NOTE: does not affect PC calculation!
%                               remove from calc using wsplice
%      plotSplicePC     0,1 indicator of whether to show splices in
%                         direction plots, default = 1
%                         NOTE: does not affect PC calculation!
%                               remove from calc using wsplice
%      scaleCurvePC     0,1 indicator of whether to scale the PC curves in
%                         dir to have maximum of +/-1
%      dirstring        cell array of strings used to label the directions,
%                         useful if specifying special directions through
%                         dir, e.g. dir = DWD dir,
%                         default = 'splice dir #' with # = 1,..,npc
%
%  parameters for violin plots
%      minN             minimum number of samples in either cluster for a
%                         splice to be plotted, default = 5
%      minP             minimum proportion of samples in either cluster for a
%                         splice to be plotted, default = 5
%
%      VPidx            length <= 2 vector of which labels to use to make
%                         violin plots,
%                         default is [2,3] when labels = []
%                         default is [1,2] when labels ~= []
%      iVPlot           0,1 indicator of whether to include violin plots,
%                         default = 1
%      iVPuse0          0,1 indicator of whether to use 0's when creating
%                         violin plots (0 counts used if iVPlot = 1),
%                         default = 1
%
%
% Output:
%   SplicePlot  - SigFuge curve plot with necessary information

%initialize default parameters for figure  
mColors = [.3 .3 .3;
           1. .3 .3; ...
           .3 .3 1.];
labels = [];
savestr = genename;
titlestr = [genename ' pile-up curves and splicing data'];

iLog = 1;
shift = 1;
iAvg = 1;
iAll = 1;
avgtype = 0;
avgP = 50;
iexbds = 1;

dir = [];
npc = 4;
minloading = 0.1;
wsplice = 0.5;
plotCurvePC = 1;
plotSplicePC = 1;
scaleCurvePC = 1;
dirstring = 'splice dir';

minN = 5;
minP = 0.05;
% VPidx;
% iVPlot;
% iVPuse0;

%mColors,labels,savestr,titlestr
%iLog,shift,iAvg,iAll,avgtype,avgP
%dir,npc,minloading,wsplice,plotCurvePC,plotSplicePC,dirstring
%minN,minP,VPidx,iVPlot,iVPuse0



if nargin > 1;   %  then paramstruct is an argument
    
    if isfield(paramstruct, 'mColors') && ~isempty(paramstruct.mColors);
        mColors = paramstruct.mColors;
    end;

    if isfield(paramstruct, 'labels');
        labels = paramstruct.labels;
    end ;

    if isfield(paramstruct, 'savestr');
        savestr = paramstruct.savestr; 
    end;

    if isfield(paramstruct, 'titlestr');
        titlestr = paramstruct.titlestr; 
    end;

    
    if isfield(paramstruct, 'iLog');
        iLog = paramstruct.iLog; 
    end ;

    if isfield(paramstruct, 'shift');
        shift = paramstruct.shift;
    end;

    if isfield(paramstruct, 'iAvg');
        iAvg = paramstruct.iAvg;
    end ;

    if isfield(paramstruct, 'iAll');
        iAll = paramstruct.iAll;
    end;
  
    if isfield(paramstruct, 'avgtype');
        avgtype = paramstruct.avgtype;
        if avgtype ~= 0;
            avgP = 5;
        end;
    end;

    if isfield(paramstruct, 'iexbds');
        iexbds = paramstruct.iexbds;
    end;

    
    if isfield(paramstruct, 'avgP');
        avgP = paramstruct.avgP;
    end;

    
    
    if isfield(paramstruct, 'dir');
        dir = paramstruct.dir; 
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

    if isfield(paramstruct, 'wsplice');
        wsplice = paramstruct.wsplice; 
    end;

    if isfield(paramstruct,'plotCurvePC');
        plotCurvePC = paramstruct.plotCurvePC; 
    end;

    if isfield(paramstruct, 'plotSplicePC');
        plotSplicePC = paramstruct.plotSplicePC; 
    end;

    if isfield(paramstruct, 'scaleCurvePC');
        scaleCurvePC = paramstruct.scaleCurvePC; 
    end;

    if isfield(paramstruct, 'dirstring');
        dirstring = paramstruct.dirstring; 
    end;

    
    
    if isfield(paramstruct, 'minN');
        minN = paramstruct.minN;
    end;

    if isfield(paramstruct, 'minP');
        minP = paramstruct.minP;
    end;


end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%load in pile-up data
depth = textread(depthpath); %#ok
depth = depth';
[d,n] = size(depth);

%labels must be nx1
%use SigFuge labels on depth if labels not provided
if length(labels) ~= n && n > 3; 
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
    disp('!!  labels not specified as nx1 vector,    !!');
    disp('!!  using default SigFuge labels           !!');
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
    labels = SigFugeLabelsPK(depth);
end;
K = max(labels);


% check whether dimensions of labels and mColors match 
if size(mColors,1) ~= K ;
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
    disp('!!  mColors must be Kx3 color matrix         !!');
    disp('!!  where K is the number of unique values   !!');
    disp('!!  in labels label vector, using default    !!');
    disp('!!  equally spaced colors on HSV colormap    !!');
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
    cmap = colormap('HSV');
    mColors = cmap(floor(linspace(1, size(cmap,1), K+1)),:);
    mColors = mColors(1:K,:);

elseif  ~all(unique(labels)==(1:size(mColors,1)));
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
    disp('!!  label vector must consist of values          !!');
    disp('!!  between 1..K (K is number of rows in mColors !!');
    disp('!!  exiting function without any output          !!');
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
    return;
end;

medColors = max(mColors-.3, 0);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load in exon data
annots = findinRefSeq({genename}, exonpath);
nexons = length(annots.start);


%load in splice site data 
splices = findinSplices(genename, splicepath, n);


%subset the set of splices, for specifics go to SpliceFilterPK.m
splices = SpliceFilterPK(splices, annots);



nsplices = length(splices);
valmat = vertcat(splices.vals);
logvalmat = log10(valmat+shift) - log10(shift);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%my guess of the genomic coordinate Chris includes in his pile-up
%curves -- not sure why some bases are missing
gcoord = annots.start(1) ;
for p = 1:nexons ;
    gcoord = [gcoord (annots.start(p)+1):annots.stop(p)] ; %#ok
end ;


%the exon boundaries for the gene
exons = annots.stop - annots.start ;
exons = [0 cumsum(exons)+1] ;
exonn = zeros(1,d) ;
for p = 1:nexons ;
  if mod(p,2) ;
    exonn((exons(p)+1):exons(p+1)) = 1 ;
  end ;
end ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%add 'introns' between genes to help with visualizing splices

ilength = ceil(d/100)*10;

%initialize with p = 1 case
exoni =   [ones((exons(2)-1)-exons(1), 1); ...
           zeros(ilength, 1)] ;
depthi =  [depth( (exons(1)+1):(exons(2)-1), : ); ...
           zeros(ilength, n)];
gcoordi = min(gcoord(exons(p)+1): ...
                (gcoord(exons(p+1))+ceil(ilength/2)), ...
              annots.start(2)-1);

          
for p = 2:(nexons-1);
    
    exoni = [exoni; ...
             ones(exons(p+1)-exons(p), 1) ; ...
             zeros(ilength, 1)]; %#ok

    depthi = [depthi; ...
              depth((exons(p)+1):exons(p+1), :); ...
              zeros(ilength, n)]; %#ok
          
      tGC = (gcoord(exons(p)+1)-floor(ilength/2)): ...
                (gcoord(exons(p+1))+ceil(ilength/2));        
      tGC = max(tGC, annots.stop(p-1)+1);
      tGC = min(tGC, annots.start(p+1)-1);
    gcoordi = [gcoordi tGC]; %#ok
    
end;


%p = nExons case
exoni =  [exoni; ...
          ones(exons(nexons+1)-(exons(nexons)-1), 1)];
depthi = [depthi; ...
          depth( exons(nexons):exons(nexons+1), : )];
  tGC = (gcoord(exons(nexons)+1)-floor(ilength/2)): ...
            gcoord(exons(nexons+1));        
  tGC = max(tGC, annots.stop(nexons-1)+1);
gcoordi = [gcoordi tGC];


%need to flip if gene is on '-' strand
if annots.dir;
    exons = exons(end) - fliplr(exons);
    depthi = flipud(depthi);
    gcoordi = fliplr(gcoordi);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%if log scale is specified, plot depth and calculate average on log-scale

logdepthi = log10(depthi+shift) - log10(shift);
if iLog;
    depthi = logdepthi;
end;

  %  use switch cases to define anonymous average function
switch avgtype
    case 0  %  median
        avgFn = @(x) quantile(x, avgP/100, 2);
    case 1  %  mean
        avgFn = @(x) mean(x, 2);
    case 2  %  winsor mean
        avgFn = @(x) mean(winsorising(x, 100-2*avgP), 2);
    case 3  %  trimmed mean
        avgFn = @(x) trimmean(x, 2*avgP, 'round', 2);
end;

%  calculate average for all curves
avg = zeros(length(depthi), K);
for k = 1:K;
  avg(:, k) = avgFn(depthi(:, labels==k));
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%calculate normalized value matricies (i.e. SS = 1)

nlogdepthi = logdepthi / ...
                sqrt( sum( n*std(logdepthi, 1, 2).^2 ) );
nlogvalmat = logvalmat / ...
                sqrt( sum( n*std(logvalmat, 1, 2).^2 ) );

fullmat = [sqrt(1-wsplice) * nlogdepthi; ...
           sqrt(wsplice) * nlogvalmat];

dplot = size(depthi, 1);

  %  defining arrows for splicing information to plot  
if ~isempty(dir);
  ddir = size(dir, 1);
  ndir = size(dir, 2);
  if ddir ~= (d+nsplices);
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
    disp('!!  specified direction vector/matrix    !!') ;
    disp('!!  dimensions differ from data!         !!') ;
    disp('!!  dir matrix must be d+dsplice x ndir. !!') ;
    disp(['!!  d+dsplice = ' num2str(d+nsplices) '               !!']) ;
    disp(['!!  ddir = ' num2str(ddir) '                            !!']) ;
    disp('!!  using default of 4 PC directions     !!') ;
    disp('!!  with curve:splices sum sq. ratio 1:1 !!') ;
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
    dir = [] ;
  else
    npc = ndir ;
    tdir = zeros(dplot,ndir) ;
    tdir(logical(exoni),:) = dir(1:d,:) ;
    if annots.dir ; 
      tdir = flipud(tdir) ;
    end ;
    dir = [tdir; dir((d+1):(d+nsplices),:)] ;
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
%plotting axes and dimensions

vaxh = axisSM(1:size(depthi, 1));
vaxv = axisSM(depthi); 
mmax = max(reshape(depthi, 1, []));
exonh1  = -.06*mmax;
exonh2  = -.09*mmax;
elabelh  = -.10*mmax;

%change plotting dimensions to account for direction plotting
vaxv(1) = -1-2*(npc+1);
PCzeros = -1-2*(1:npc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%actually plotting the figure

figure(1);
clf;
hold on;
ax = [vaxh(1:2) vaxv(1:2)];
axis(ax);

%add exon boundaries
icolor = hex2dec({'5B' 'A4' 'CB'})' / 255;
ecolor = hex2dec({'00' '40' '62'})' / 255;

rectangle('Position', [0, -.08*mmax, ...
                      (nexons-1)*ilength+exons(nexons+1), 0.01*mmax], ...
          'FaceColor', icolor, ...
          'EdgeColor', icolor);
for p = 1:nexons;
  rectangle('Position', [exons(p)+(p-1)*ilength, -.10*mmax, ...
                        exons(p+1)-exons(p), .05*mmax], ...
            'FaceColor', ecolor, ...
            'EdgeColor', ecolor);
end;
text(-d/150, (exonh1+exonh2)/2, ...
     annots.chr, ...
     'HorizontalAlignment', 'right', ...
     'VerticalAlignment', 'middle', ...
     'FontSize', 10);
    

%background rectangle for L4, L6 and exon boundaries
Xeps = (ax(2)-ax(1))*.02;
for i = 1:npc;
  rectangle('Position', [ax(1)+Xeps,-2*(i+1),(ax(2)-ax(1))-2*Xeps,2], ...
            'FaceColor', [.9 .9 .9]-mod(i,2)*[.2 .2 .2], ...
            'EdgeColor', [.9 .9 .9]-mod(i,2)*[.2 .2 .2], ...
            'curvature', 0);
end;

   
%adding genomic coordinates
if iexbds;
    for p = 1:nexons ;
        if annots.dir ;
            pStart = annots.stop(nexons-p+1) ;
            pEnd = annots.start(nexons-p+1) ;
        else
            pStart = annots.start(p) ;
            pEnd = annots.stop(p) ;
        end ;
        text(exons(p)+(p-1)*ilength, ...
            elabelh, ...
            num2str(pStart), ...
            'rotation',60, ...
            'horizontalalignment','right', ...
            'verticalalignment','top', ...
            'FontSize',9) ;      
        text(exons(p+1)+(p-1)*ilength-20, ...
            elabelh, ...
            num2str(pEnd), ...
            'rotation',60, ...
            'horizontalalignment','right', ...
            'verticalalignment','top', ...
            'FontSize',9) ;  
    end ;
end;

  %  bottom lines for PC loadings plots
line(vec2matSM([0; length(depthi)], npc), vec2matSM(PCzeros, 2), ...
     'color', 'k', ...
     'linewidth', 2) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Calculating the splice plots for PC directions %%%%%%
   
  %  computing splice arrow start and end positions 
  %  p1,p2 = position1 (start) --> position2 (end)
if annots.dir ;
  p1 = arrayfun(@(x) sum(x<gcoordi),[splices.stop]) ;
  p2 = arrayfun(@(x) sum(x<gcoordi),[splices.start]) ;
else
  p1 = arrayfun(@(x) sum(x>gcoordi),[splices.start]) ;
  p2 = arrayfun(@(x) sum(x>gcoordi),[splices.stop]) ;        
end ;
p1(p1==0) = -ilength/2 ;
  %  extend some splices before gene model
p2(p2==length(gcoordi)) = length(gcoordi) + ilength/2 ;
  %  extend some splice beyond gene model
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%actually plotting all PC or user specified directions

for i = 1:npc ;

    ic = (1-mod(i,2))*[.5 0 .5] + mod(i,2)*[0 .5 0] ;
    axis(axis) ;
    if plotSplicePC && max(abs(mdirval((dplot+1):end,i)))>0 ;
        PCarrowH = mdirval((dplot+1):end,i) + PCzeros(i) ;
        vp1 = [p1' PCarrowH] ; vp2 = [p2' PCarrowH] ;
        vp1 = vp1((abs(mdirval((dplot+1):end,i))>=minloading),:) ;
        vp2 = vp2((abs(mdirval((dplot+1):end,i))>=minloading),:) ;
        if ~isempty(vp1) ;
            PCarrows = arrow(vp1,vp2,'Length',4,'TipAngle',30,'Width',1, ...
                          'EdgeColor',ic,'FaceColor',ic) ; %#ok
        end ;
    end ;

      %plot the curve loadings (assuming, again, [curves; splices] format
    if plotCurvePC && max(abs(mdirval(1:dplot,i)))>0 ;
        maxPCload = max(abs(mdirval(1:dplot, i)));
        if scaleCurvePC || maxPCload > 1;
            mdirval(1:dplot, i) = mdirval(1:dplot, i) / ...
                                          maxPCload;
        end;
        PCcurveH = mdirval(1:dplot, i) + PCzeros(i) ;
        plot(PCcurveH,'-','color',ic) ;
    end ;

end ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%splice plots for top plot

SpliceCnts = zeros(size(valmat,1), K);
for k = 1:K;
  SpliceCnts(:,k) = sum(valmat(:,labels==k)>0, 2);
end;

cn = zeros(1, K);
for k = 1:K;
  cn(k) = sum(labels==k);
end;
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

fInterest = cellfun(@(x) any(x>=min(floor(cn*minP), minN)), ...
                    num2cell(SpliceCnts, 2));
plotACT = splices(fInterest);

plotvalmat = vertcat(plotACT.vals);
plotvalmat = log10(plotvalmat+shift) - log10(shift);

mArrowH = zeros(size(plotvalmat,1), K);
for k = 1:K;
    mArrowH(:,k) = avgFn(plotvalmat(:,labels==k));
end ;

%add all pile-up curves to top panel
if iAll;
    for k = 1:K;
        plot(depthi(:,labels==k), '-', ...
             'color', mColors(k,:));
    end;
end;

%add average pile-up curves to top panel
if iAvg;
    for k = 1:K;
        plot(avg(:,k), '-', ...
             'color', medColors(k,:), ...
             'linewidth', 2);
    end;
end;



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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%labeling X,Y axes and adding annotations before printing

%removing X-axis
set(gca, 'XTick', []);

%re-labeling Y-axis
set(gca, 'YTick', [sort(PCzeros) 0:floor(vaxv(2))]);
temp = get(gca, 'YTickLabel');
temp = mat2cell(temp, ones(size(temp,1),1), size(temp,2));
if isempty(dir);
    dirtype = 'splice PC';
else
    dirtype = dirstring;
end;

if length(dirstring) == npc;
    temp(1:npc) = fliplr(dirstring);
else 
    for i = 1:npc;
        temp(npc+1-i) = {[dirtype num2str(i)]};
    end;
end;
set(gca, 'YTickLabel', temp);


%adding X,Y axis labels and title to plot
xlabel('genomic position');
ylabel('log_{10} (counts+1)');
if ~isempty(dir);
    title(titlestr);
else
    title([genename ' pile-up curves and splice PC loadings']);
end;

%printing cluster/group sizes to upper right of figure
ax = axis;
for k = 1:K;
    text(ax(1)+(ax(2)-ax(1))*.95, ...
         ax(3)+(ax(4)-ax(3))*(1.06-k*.03), ...
         ['cluster ' num2str(k) ' size = ' num2str(cn(k))], ...
         'color',mColors(k,:), ...
         'FontSize', 10);
end;


orient landscape ;
print('-dpdf',savestr) ;   



