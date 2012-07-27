% afilename='/usr/scratch2/ulrich/DrosophilaEmbryo/NuclearDynamics/hisB.tif';
% outfolder='/usr/scratch2/ulrich/DrosophilaEmbryo/NuclearDynamics/';
afilename='/usr/scratch2/ulrich/DrosophilaEmbryo/Wavy/Wavy.tif';
datafolder='/usr/scratch2/ulrich/DrosophilaEmbryo/Wavy/';
outfolder='/usr/scratch2/ulrich/DrosophilaEmbryo/Wavy/fine/';

picnumMax=164;
r_Zr = [24 28 32 36 41 50 75 100 140 200 400 800];
r_Zr = [20:1:100 104:4:600];
colorRange=128;
colors=jet(colorRange);colors(colorRange+1:2000,:)=colors(colorRange*ones(2000-colorRange,1),:);

ZrNpNe=cell(picnumMax,3);
ZrAbs=zeros(picnumMax,length(r_Zr));
meanEdgLengths=zeros(picnumMax,1);

picnum=13;
%for picnum=1:picnumMax

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

%figure(2);hold on;triplot(DT);hold off;
%edgs=edges(DT);


xInside=x(InsideQ);
yInside=y(InsideQ);


vx=x(edgsInside');
vy=y(edgsInside');
%figure(2);hold on;plot(vx,vy,'r-');hold off;

CH=convexHull(DT);

edgLengths2 = ( vx(1,:)-vx(2,:) ).^2 + ( vy(1,:)-vy(2,:) ).^2;
meanEdgLength=(mean(edgLengths2.^(-2)))^(-1/4)
meanEdgLengths(picnum)=meanEdgLength;
%numberOfEntities=length(x);


figure(3);clf('reset');
PlotVoronoi;
%print(gcf,'-dpng',[outfolder 'out' num2str(picnum) '.png'],'-r200');


%figure(51);clf('reset');
%PlotZr
%print(gcf,'-dpng',[outfolder 'Zr' num2str(picnum) '.png'],'-r200');

%end
%if all(ZrAbs~=0)
%    disp(['saving Zr-Data to ' outfolder 'outZr.txt']);
%    save([outfolder 'outZr.txt'],'ZrNpNe', 'ZrAbs', 'meanEdgLengths', 'r_Zr');
%end