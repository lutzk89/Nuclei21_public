% afilename='/usr/scratch2/ulrich/DrosophilaEmbryo/NuclearDynamics/hisB.tif';
% outfolder='/usr/scratch2/ulrich/DrosophilaEmbryo/NuclearDynamics/';
afilename='/usr/scratch2/ulrich/DrosophilaEmbryo/Wavy/Wavy.tif';
datafolder='/usr/scratch2/ulrich/DrosophilaEmbryo/Wavy/';
outfolder='/usr/scratch2/ulrich/DrosophilaEmbryo/Wavy/wavePlot/';

picnumMax=164;
r_Zr = 1;
colorRange=128;
colors=jet(colorRange);colors(colorRange+1:2000,:)=colors(colorRange*ones(2000-colorRange,1),:);
numberOfSlices=15;

SlicePartilceN      = zeros(picnumMax,numberOfSlices);
SliceDensities      = zeros(picnumMax,numberOfSlices);
SliceMeanEdgLengths = zeros(picnumMax,numberOfSlices);
SliceOrderParameter = zeros(picnumMax,numberOfSlices);
SliceF6Neighbors    = zeros(picnumMax,numberOfSlices);
SliceDeansityMaxSlope = zeros(picnumMax,2);

picnum=62;
LastMaxSlope=0;
for picnum=1:picnumMax
%for picnum=50:95

imSav0=double(imread(afilename,'tif',picnum));
%imSav0=ones(800,1100);
figure(1);hold off;image(imSav0*.29);colormap bone;
imSav1=64*imSav0/max(max(imSav0));


SAVMAT=load([datafolder 'outPos' num2str(picnum) '.txt'], '-ascii');
x=SAVMAT(:,1)';
y=SAVMAT(:,2)';
n=length(x);

%createRandomConf;
figure(1)
hold on
plot(x,y,'r.')
hold off


DT=DelaunayTri(x',y');
[ InsideQ, edgsInside ]=findInsideCH(DT,40);

[V,R] = voronoiDiagram(DT);
InsideIndices=find(InsideQ);

edgs=edges(DT);
S=sparse(edgs(:,1),edgs(:,2),1,n,n);
SSym0=S+S';

NNeighbors=full(sum(SSym0,1))';
Inside6=(InsideQ & NNeighbors==6);
Inside5=(InsideQ & NNeighbors<=5);
Inside7=(InsideQ & NNeighbors>=7);

xInside=x(InsideQ);
yInside=y(InsideQ);
SSym0Inside=SSym0(InsideIndices,InsideIndices);
[Zr,Np,Ne]=calcZrS2( xInside, yInside, SSym0Inside , r_Zr);
ZrAbs=abs(Zr./Ne);

xMax=size(imSav0,2);
yMax=size(imSav0,1);


for ii=1:length(InsideIndices)
    Ind   = InsideIndices(ii);
    slice = floor(x(Ind)/xMax*numberOfSlices) + 1;
    SlicePartilceN(picnum,slice) = SlicePartilceN(picnum,slice) + 1;
    SliceDensities(picnum,slice) = SliceDensities(picnum,slice) + 1/polyarea(V(R{Ind},1),V(R{Ind},2));
    neighborInds = find(SSym0(:,Ind));
    SliceMeanEdgLengths(picnum,slice) = SliceMeanEdgLengths(picnum,slice) + mean(sqrt( (x(Ind)-x(neighborInds)).^2 + (y(Ind)-y(neighborInds)).^2 ));
    SliceOrderParameter(picnum,slice) = SliceOrderParameter(picnum,slice) + ZrAbs(ii);
    if length(neighborInds)==6
        SliceF6Neighbors(picnum,slice) = SliceF6Neighbors(picnum,slice) + 1;
    end
end

for slice=1:numberOfSlices
    SliceDensities(picnum,slice)      = SliceDensities(picnum,slice)      / SlicePartilceN(picnum,slice);
    SliceMeanEdgLengths(picnum,slice) = SliceMeanEdgLengths(picnum,slice) / SlicePartilceN(picnum,slice);
    SliceOrderParameter(picnum,slice) = SliceOrderParameter(picnum,slice) / SlicePartilceN(picnum,slice);
    SliceF6Neighbors(picnum,slice)    = SliceF6Neighbors(picnum,slice)    / SlicePartilceN(picnum,slice);
end
if picnum==50
    LastMaxSlope=0;
end
if LastMaxSlope==0
    [SliceDeansityMaxSlope(picnum,1), SliceDeansityMaxSlope(picnum,2)] = max(SliceDensities(picnum,2:numberOfSlices) - SliceDensities(picnum,1:(numberOfSlices-1)));
    LastMaxSlope=SliceDeansityMaxSlope(picnum,2);
else
    checkers=max(LastMaxSlope-1 ,2):min(LastMaxSlope+2 ,numberOfSlices);
    [SliceDeansityMaxSlope(picnum,1), SliceDeansityMaxSlope(picnum,2)] = max(SliceDensities(picnum,checkers) - SliceDensities(picnum,checkers-1));
    SliceDeansityMaxSlope(picnum,2)=SliceDeansityMaxSlope(picnum,2)+checkers(1)-2;
    LastMaxSlope=SliceDeansityMaxSlope(picnum,2);
end


figure(2);clf
plotRange=2:(numberOfSlices-1);
%%%%%%%%%%%%%%%%% Plot 1 %%%%%%%%%%%%%%%%%
subplot(2,2,1);
hold off; image(imSav0*.29); colormap bone;
hold on; 
plot(xInside,yInside,'r.');
for ii=0:xMax/numberOfSlices:xMax
    plot([ii ii],[0 yMax]);
end
hold off;


plotRangeMax=[SliceDeansityMaxSlope(picnum,2) SliceDeansityMaxSlope(picnum,2)+1];
%%%%%%%%%%%%%%%%% Plot 2 %%%%%%%%%%%%%%%%%
subplot(2,2,2);
plot(plotRange,SliceDensities(picnum,plotRange));axis([0 numberOfSlices 0 0.0015]);
hold on
plot(plotRangeMax,SliceDensities(picnum,plotRangeMax),'r');
hold off
title('Local Density')

%%%%%%%%%%%%%%%%% Plot 3 %%%%%%%%%%%%%%%%%
subplot(2,2,3);
plot(plotRange,SliceOrderParameter(picnum,plotRange));axis([0 numberOfSlices 0 1]);
if 50<=picnum<=95
    hold on
    plot(plotRangeMax,SliceOrderParameter(picnum,plotRangeMax),'r');
    hold off
end
title('Local Order Parameter')


%%%%%%%%%%%%%%%%% Plot 4 %%%%%%%%%%%%%%%%%
subplot(2,2,4);
plot(plotRange,SliceF6Neighbors(picnum,plotRange));axis([0 numberOfSlices 0 1]);
if 50<=picnum<=95
    hold on
    plot(plotRangeMax,SliceF6Neighbors(picnum,plotRangeMax),'r');
    hold off
end
title('Fraction of Nuclei with 6 NN')

% %%%%%%%%%%%%%%%%% Plot 4 %%%%%%%%%%%%%%%%%
% subplot(2,2,4);
% plot(plotRange,SliceMeanEdgLengths(picnum,plotRange));axis([0 numberOfSlices 0 60]);
% title('Mean Edge Length')

print(gcf,'-dpng',[outfolder 'WaveNNmaxSlope' num2str(picnum) '.png'],'-r200');
end
save([outfolder 'outWave.txt'], 'SliceDensities', 'SliceMeanEdgLengths', 'SliceOrderParameter', 'SliceF6Neighbors', 'SlicePartilceN','SliceDeansityMaxSlope');

plotRangeTime=50:95;
figure(3);plot(plotRangeTime,SliceDeansityMaxSlope(plotRangeTime,2));
Times = 60:85;
poss  = SliceDeansityMaxSlope(Times,2)';
timeBar = mean(Times);
posBar  = mean(poss);
vel  = sum( (Times-timeBar).*(poss-posBar) ) / sum( (Times-timeBar).^2 );
pos0 = posBar - vel*timeBar;
hold on
plot(Times,pos0+vel*Times,'r')
hold off

