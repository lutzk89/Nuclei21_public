% afilename='/usr/scratch2/ulrich/DrosophilaEmbryo/NuclearDynamics/hisB.tif';
% outfolder='/usr/scratch2/ulrich/DrosophilaEmbryo/NuclearDynamics/';
afilename='/usr/scratch2/ulrich/DrosophilaEmbryo/Wavy/Wavy.tif';
outfolder='/usr/scratch2/ulrich/DrosophilaEmbryo/Wavy/';

thresholdFAC1=1.0;
thresholdFAC2=1.2;
maskBlackFAC=0.3;
preMaskFAC=0.12;
high_limitFAC=3.5e-4;

picnum=100;
%for picnum=1:164

imSav=double(imread(afilename,'tif',picnum));
figure(1);hold off;image(imSav*.29);colormap bone;
%imSav=imSav-imSavMean;
%imSav=log(imSav+.01);
imSav=log(log(imSav+1)+1);
imSavMean=mean(mean(imSav));%*.65;
imSav=imSav-imSavMean;
preMask=circle(6);
maskWhite=circle(2);
maskBlack=circle(9);

h=bp_sample(.045,0.18,80);
im=filter2s(h,imSav);
figure(2);imagesc(im);colorbar;
threshold=thresholdFAC1*sqrt(mean(mean(im.^2)) - mean(mean(im))^2);
numthreshold=5;

[x,y,Num,Sum]=entities2(im,threshold,numthreshold,preMask,maskWhite,maskBlack);

%figure(2);hold on;plot(x,y,'b.');hold off;
DT=DelaunayTri(x',y');
%figure(2);hold on;triplot(DT);hold off;
edgs=edges(DT);

n=length(x);S=sparse(edgs(:,1),edgs(:,2),1,n,n);
SSym=S+S'+sparse(1:n,1:n,1);
%boundaryIndices=convexHull(DT);
boundaryIndices=zeros(n,1);boundaryIndices(convexHull(DT))=1;
boundaryIndices2=SSym*boundaryIndices;
InsideIndices=~boundaryIndices2;

xinside=x(InsideIndices);
yinside=y(InsideIndices);
figure(2);hold on;plot(xinside,yinside,'b.');hold off;

boundaryIndices2s=find(boundaryIndices2~=0);
removeEdgs=(edgs==boundaryIndices2s(1));
for i=2:length(boundaryIndices2s)
    removeEdgs=removeEdgs+(edgs==boundaryIndices2s(i));
end
removeEdgs1=(removeEdgs(:,1)+removeEdgs(:,2)~=0);
edgsInside=edgs;
edgsInside(removeEdgs1,:)=[];
vx=x(edgsInside');
vy=y(edgsInside');
figure(2);hold on;plot(vx,vy,'r-');hold off;


edgLengths2 = ( vx(1,:)-vx(2,:) ).^2 + ( vy(1,:)-vy(2,:) ).^2;
meanEdgLength=(mean(edgLengths2.^(-2)))^(-1/4)
numberOfEntities=length(x);




low_limit=.03;
%high_limit=cutOff((16./meanEdgLength)^2, 0.075, 1);
high_limit=cutOff(high_limitFAC*numberOfEntities, 0.075, 1);

h=bp_sample(low_limit,high_limit,100);%-0.4*bp_sample(0.0,.01,100);
%h2=bp_sample(0,.4,25);

im=filter2s(h,imSav);
%im=cutoff(im,0.,2*sqrt(mean(mean(im.^2)) - mean(mean(im))^2));
figure(5);imagesc(im);colorbar
%figure(5);imagesc(cutOff(im,0,1000.0));colorbar;
%ims=im(YT(channel+1):YB(channel+1),XL(channel+1):XR(channel+1));

%imGrob=filter2(bp_sample(0.0,.01,60),imSav);
%figure(5);imagesc(imGrob);colorbar
%im=im-imGrob;
%figure(5);imagesc(im);colorbar

%circleSize=round(meanEdgLength*.2);
%maskBlack=circle(circleSize+6);
%maskWhite=circle(circleSize-5);

maskBlack=circle(meanEdgLength*maskBlackFAC);
%maskWhite=circle(round(meanEdgLength*.04));
maskWhite=[1 1; 1 1];
%preMask=circle(meanEdgLength*.1);
preMask=circle(meanEdgLength*preMaskFAC);
%threshold=42000/meanEdgLength^2%;25;%imSavMean*.3%20;%15-imSavMean;
%threshold=thresholdFAC2* (sqrt(mean(mean(im.^2)) - mean(mean(im))^2));%/(meanEdgLength^.2);%6.5/(meanEdgLength^.4);
threshold=2.1*thresholdFAC2*(mean(mean(cutOff(im,0,1000.0))));
numthreshold=1./high_limit;%meanEdgLength*1.2;
[x,y,Num,Sum]=entities2(im,threshold,numthreshold,preMask,maskWhite,maskBlack);


figure(1)
hold on
plot(x,y,'r.')
hold off


%saveas(gcf,[outfolder '/out' num2str(picnum) '.png']);
%end