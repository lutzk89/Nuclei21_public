afilename='/usr/scratch2/ulrich/DrosophilaEmbryo/NuclearDynamics/hisB.tif';
outfolder='/usr/scratch2/ulrich/DrosophilaEmbryo/NuclearDynamics/';


picnum=2;
for picnum=1:120

imSav=double(imread(afilename,'tif',picnum));
figure(1);hold off;image(imSav*.29);colormap bone;
imSavMean=mean(mean(imSav))*.75;
imSav=imSav-imSavMean;
maskWhite=circle(2);
maskBlack=circle(9);

h=bp_sample(.03,0.2,45);
im=filter2(h,imSav);
figure(2);imagesc(im);colorbar;
threshold=28;
numthreshold=5;

[x,y,Num,Sum]=entities2(im,threshold,numthreshold,circle(7),maskWhite,maskBlack);

figure(2);hold on;plot(x,y,'b.');hold off;
edgs=edges(DelaunayTri(x',y'));
edgLengths = sqrt( ( x(edgs(:,1))-x(edgs(:,2)) ).^2 + ( y(edgs(:,1))-y(edgs(:,2)) ).^2 );
meanEdgLength=mean(edgLengths)




low_limit=.03;
high_limit=(32.2/meanEdgLength)^3.5;
h=bp_sample(low_limit,high_limit,60)-2*bp_sample(0.0,.01,60);
%h2=bp_sample(0,.4,25);

im=filter2(h,imSav);
figure(4);imagesc(im);colorbar
%ims=im(YT(channel+1):YB(channel+1),XL(channel+1):XR(channel+1));

%imGrob=filter2(bp_sample(0.0,.01,60),imSav);
%figure(5);imagesc(imGrob);colorbar
%im=im-imGrob;
%figure(5);imagesc(im);colorbar

circleSize=round(meanEdgLength*.2);
%maskBlack=circle(circleSize+6);
%maskWhite=circle(circleSize-5);

maskBlack=circle(round(meanEdgLength*.3));
maskWhite=circle(round(meanEdgLength*.04));
preMask=circle(round(meanEdgLength*.13));
%threshold=42000/meanEdgLength^2%;25;%imSavMean*.3%20;%15-imSavMean;
threshold=sqrt(mean(mean(im.^2)) - mean(mean(im))^2)*3.7/(meanEdgLength^.25);%6.5/(meanEdgLength^.4);
numthreshold=1/high_limit;%meanEdgLength*1.2;
[x,y,Num,Sum]=entities2(im,threshold,numthreshold,preMask,maskWhite,maskBlack);


figure(1)
hold on
plot(x,y,'r.')
hold off


saveas(gcf,[outfolder '/out' num2str(picnum) '.png']);
end