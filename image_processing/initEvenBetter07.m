% afilename='/usr/scratch2/ulrich/DrosophilaEmbryo/NuclearDynamics/hisB.tif';
% outfolder='/usr/scratch2/ulrich/DrosophilaEmbryo/NuclearDynamics/';
% afilename='/usr/scratch2/ulrich/DrosophilaEmbryo/Wavy/Wavy.tif';
% outfolder='/usr/scratch2/ulrich/DrosophilaEmbryo/Wavy/';
afilename='/usr/scratch2/ulrich/DrosophilaEmbryo/better_data/d07.tif';
outfolder='/usr/scratch2/ulrich/DrosophilaEmbryo/better_data/d07/';

fileInfo=imfinfo(afilename,'tif');

%[XXX YYY]=ndgrid(1:fileInfo(1).Height,1:fileInfo(1).Width); repaintIndices=find((YYY-835)*90>XXX*415); clear XXX YYY;


colorRange=128;
colors=jet(colorRange);colors(colorRange+1:2000,:)=colors(colorRange*ones(2000-colorRange,1),:);

thresholdSkew1 = 1.2;
thresholdSkew2 = 1.0;
preMaskFAC         = 0.1;
maskBlackFAC       = 0.1;
maskBlackRadiusMax = 18;
high_limitFAC = 3.;
%high_limitFAC=3.5e-4;


%picnumMax = 322;
%picnumMax = 176;
picnumMax=length(fileInfo);

picnum = 40;
%for picnum=50:picnumMax
for picnum=1:picnumMax

imSav0=double(imread(afilename,'tif',picnum));
%imSav0(repaintIndices)=78.;

figure(1);hold off;image(imSav0*.29);colormap bone;
set(1,'Position',[10 620 size(imSav0,2)/2. size(imSav0,1)/2.])

imSav1=64*imSav0/max(max(imSav0));
%imSav=imSav-imSavMean;
%imSav=log(imSav+.01);
imSav=log(log(imSav0+1)+1);
imSavMean=mean(mean(imSav));%*.65;
imSav=imSav-imSavMean;


%h=bp_sample(.031,0.18,90);
%h=bp_sample(.031,0.09,90);
h=bp_sample(.031,0.06,90);
im=filter2s(h,imSav);
figure(2);imagesc(im);colorbar;
set(2,'Position',[610 620 size(imSav0,2)/2. size(imSav0,1)/2.])

circleSize100filter = 25;
filterMaxNorm
%threshold    = 0.01;
threshold    = mean(im100(im100>0));
numthreshold = 2;
preMask   = circle(10);
maskBlack = circle(3);
maskWhite = circle(2);

[x,y,Num,Sum]=entities2(im100,threshold,numthreshold,preMask,maskWhite,maskBlack);


%threshold=thresholdFAC1*sqrt(mean(mean(im.^2)) - mean(mean(im))^2);
%threshold=thresh(im, thresholdSkew1);
%numthreshold=5;


%[x,y,Num,Sum]=entities2(im,threshold,numthreshold,preMask,maskWhite,maskBlack);

%figure(2);hold on;plot(x,y,'b.');hold off;
DT=DelaunayTri(x',y');
%figure(2);hold on;triplot(DT);hold off;
edgs=edges(DT);
[ InsideQ, edgsInside ]=findInsideCH(DT,40);

xinside=x(InsideQ);
yinside=y(InsideQ);
figure(2);hold on;plot(xinside,yinside,'b.');hold off;

vx=x(edgsInside');
vy=y(edgsInside');
%figure(2);hold on;plot(vx,vy,'r-');hold off;

CH=convexHull(DT);
figure(2);hold on;plot(x(CH), y(CH), 'b');hold off;

edgLengths2 = ( vx(1,:)-vx(2,:) ).^2 + ( vy(1,:)-vy(2,:) ).^2;
meanEdgLength=(mean(edgLengths2.^(-1)))^(-1/2)
numberOfEntities=length(x);


% figure(3);clf('reset');
% 
% [V,R] = voronoiDiagram(DT);
% InsideIndices=find(InsideQ);
% areas=zeros(length(x),1);
% for i=1:length(InsideIndices)
%     Ind=InsideIndices(i);
%     areas(Ind)=polyarea(V(R{Ind},1),V(R{Ind},2));
% end
% 
% %MaxArea=max(areas);
% MaxArea=meanEdgLength^2/4*pi;
% areas=cutOff(areas/MaxArea,0,2);
% 
% 
% subplot(1,2,1);image(imSav1);
% axis([0 size(imSav0,2) 0 size(imSav0,1)]);axis equal;axis image;
% subplot(1,2,2);image(imSav1);
% axis([0 size(imSav0,2) 0 size(imSav0,1)]);axis equal;axis image;
% hold on;
% for i=1:length(InsideIndices)
%     Ind=InsideIndices(i);
%     patch(V(R{Ind},1),V(R{Ind},2),64*areas(Ind));
% end
% hold off;
% 





%low_limit=.03;
low_limit=.03;
high_limit=cutOff((high_limitFAC/meanEdgLength), 0.045, 1);
%high_limit=cutOff((high_limitFAC/meanEdgLength)^2, 0.075, 1);
%high_limit=cutOff(high_limitFAC*numberOfEntities, 0.075, 1);

h=bp_sample(low_limit,high_limit,100);%-0.4*bp_sample(0.0,.01,100);
%h2=bp_sample(0,.4,25);

im=filter2s(h,imSav);
%im=cutoff(im,0.,2*sqrt(mean(mean(im.^2)) - mean(mean(im))^2));
figure(5);imagesc(im);colorbar
set(5,'Position',[610 420 size(imSav0,2)/2. size(imSav0,1)/2.])

circleSize100filter = meanEdgLength*.6;
filterMaxNorm
threshold    = mean(im100(im100>0));%/10;
numthreshold = 2;
preMask   = circle(meanEdgLength*preMaskFAC);
maskBlackRadius = min(meanEdgLength*maskBlackFAC, maskBlackRadiusMax);
maskBlack=circle(maskBlackRadius);
maskWhite = circle(2);

[x,y,Num,Sum]=entities2(im100,threshold,numthreshold,preMask,maskWhite,maskBlack);


%figure(5);imagesc(cutOff(im,0,1000.0));colorbar;
%ims=im(YT(channel+1):YB(channel+1),XL(channel+1):XR(channel+1));

%imGrob=filter2(bp_sample(0.0,.01,60),imSav);
%figure(5);imagesc(imGrob);colorbar
%im=im-imGrob;
%figure(5);imagesc(im);colorbar

%circleSize=round(meanEdgLength*.2);
%maskBlack=circle(circleSize+6);
%maskWhite=circle(circleSize-5);




% 
% maskBlackRadius = min(meanEdgLength*maskBlackFAC, maskBlackRadiusMax);
% maskBlack=circle(maskBlackRadius);
% %maskWhite=circle(round(meanEdgLength*.04));
% maskWhite=[1 1; 1 1];
% %maskWhite=[1 1 1; 1 1 1; 1 1 1];
% %preMask=circle(meanEdgLength*.1);
% preMask=circle(meanEdgLength*preMaskFAC);
% %threshold=42000/meanEdgLength^2%;25;%imSavMean*.3%20;%15-imSavMean;
% %threshold=thresholdFAC2* (sqrt(mean(mean(im.^2)) - mean(mean(im))^2));%/(meanEdgLength^.2);%6.5/(meanEdgLength^.4);
% %threshold=2.1*thresholdFAC2*(mean(mean(cutOff(im,0,1000.0))));
% threshold=thresh(im,thresholdSkew2 + max((meanEdgLength-50)/40, 0));
% numthreshold=2./high_limit;%meanEdgLength*1.2;
% [x,y,Num,Sum]=entities2(im,threshold,numthreshold,preMask,maskWhite,maskBlack);
% 





figure(1)
hold on
plot(x,y,'r.');
hold off
set(gcf, 'PaperPositionMode', 'auto');
print(gcf,'-dpng',[outfolder 'pic' num2str(picnum) '.png'],'-r200')

SAVMAT=[x' y' Num' Sum']; save([outfolder 'outPos' num2str(picnum) '.txt'],'SAVMAT','-ascii', '-tabs');

% 
% DT=DelaunayTri(x',y');
% %figure(2);hold on;triplot(DT);hold off;
% edgs=edges(DT);
% [ InsideQ, edgsInside ]=findInsideCH(DT,40);
% 
% xInside=x(InsideQ);
% yInside=y(InsideQ);
% 
% 
% vx=x(edgsInside');
% vy=y(edgsInside');
% %figure(2);hold on;plot(vx,vy,'r-');hold off;
% 
% CH=convexHull(DT);
% 
% %edgLengths2 = ( vx(1,:)-vx(2,:) ).^2 + ( vy(1,:)-vy(2,:) ).^2;
% %meanEdgLength=(mean(edgLengths2.^(-2)))^(-1/4)
% %numberOfEntities=length(x);
% 
% 
% figure(3);clf('reset');
% PlotVoronoi;
% %print(gcf,'-dpng',[outfolder 'out' num2str(picnum) '.png'],'-r200');

end