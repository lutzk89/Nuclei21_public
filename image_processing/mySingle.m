FigNums = [1 2 3 100 5 6 100];
%close(FigNums);
filename='data/Lutz/area3.tif';
colorRange=128;
colors=jet(colorRange);colors(colorRange+1:2000,:)=colors(colorRange*ones(2000-colorRange,1),:);


%% Parameters for the first search 

thresholdSkew1  = 0.9;
thresholdFAC1   = 0.9;
low_limit1      = .001;

FFTSpectrumTilt = 120.+0*4.;
% If the wrong "image wavelength" is found in figure 2, this parameter tilts (skews) the
% Fourier specturm. 
% Larger values favors larger wavelengths

FFTLocalMaxKernel = 1;
%FFTLocalMaxKernel = [.25,.5,.25];


%% Parameters for the second search


thresholdFAC2  = .9;
% Threshold for the image binarization
% smaller value: find more nuclei,  connect closeby entities
%  larger value: find fewer nuclei, separate closeby entities

                        

high_limitFAC0 = 4.5;
%high_limitFAC  = 3.;
high_limitFAC  = 1.5;
% larger value: less smoothing, closeby entities are better separated                         
% smaller value: more smoothing, closeby entities are more often 
% considered as one entity


%low_limit2     = .015;
low_limit2     = .35;  
% larger  value: "flatten out" inhomogeneities in illumination
% smaller value: leave image more the way it is


preMaskFAC     = 0.25;   
maskBlackRadiusMax = 25;3
maskBlackFAC   = 0.1;   
%  larger value: closeby entities are better separated
% smaller value: closeby entities are more often 
% considered as one entity
                        
                        
thresholdSkew2 = 0.0;   % not very important                        
subtract100fac = 0.6;   % not very important
LargeNucDist   = 10.0;


%% Read Image 
imSav0=imread(filename,'tif');
figure, imshow(imSav0)
%if FigNums(1)>0
%    figure(FigNums(1));hold off;image(imSav0*.29);colormap bone;
%    set(FigNums(1),'Position',[10 630 size(imSav0,2)/2. size(imSav0,1)/2.])
%end

imSav1 = 64*imSav0/max(max(imSav0));
imSav=imSav0-imSavMean;
%imSav=log(imSav+.01);
%imSav     = log(log(imSav0+1)+1);
%imSav     = log(imSav0+1);
imSavMean = mean(mean(imSav));%*.65;
%imSav     = imSav-imSavMean;
%figure(41);imagesc(imSav);colorbar;

%% Find prevailing wavelength from Fourier Spectrum

meanEdgLength0 = findFrqLog(imSav, 50, FFTSpectrumTilt, FFTLocalMaxKernel, 1., FigNums(2));


%% Filter image for first search
disp(['prevailing image wave length = ' num2str(meanEdgLength0)]);
high_limit0 = cutOff((high_limitFAC0/meanEdgLength0), 0.032, 1);
h           = bp_sample(low_limit1, high_limit0, 450);
%h           = bp_sample(.05, high_limit0, 450);

im = filter2fft(h,imSav);

% if FigNums(3)>0
%     figure(FigNums(3));imagesc(im);colorbar;
%     set(FigNums(3),'Position',[610 630 size(imSav0,2)/2. size(imSav0,1)/2.])
% end




%% Binarize image and find entities
circleSize100filter = 25;
thresholdSkew = thresholdSkew1;
filterMaxNorm
%threshold    = 0.01;
figure, imshow(im100)
mine=im2bw(im100,graythresh(im100));% if FigNums(6)>0
%     figure(FigNums(6));subplot(2,2,3);imagesc(im100);colorbar
% end

figure, imshow(mine)
threshold    = mean(im100(im100>0)) * thresholdFAC1;
numthreshold = 2;
preMask   = circle(3);
maskBlack = circle(3);
maskWhite = circle(2);


[x,y,Num,Sum]=entities2(im100,threshold,numthreshold,preMask,maskWhite,maskBlack,[FigNums(4), -1]);


%threshold=thresholdFAC1*sqrt(mean(mean(im.^2)) - mean(mean(im))^2);
%threshold=thresh(im, thresholdSkew1);
%numthreshold=5;


%[x,y,Num,Sum]=entities2(im,threshold,numthreshold,preMask,maskWhite,maskBlack);


DT=DelaunayTri(x',y');
%if FigNums(3)>0
%    figure(FigNums(3));hold on;plot(x,y,'b.');hold off;
%    figure(FigNums(3));hold on;triplot(DT);hold off;
%end
edgs = edges(DT);
[ InsideQ, edgsInside ] = findInsideCH(DT,40);

xinside=x(InsideQ);
yinside=y(InsideQ);
% if FigNums(3)>0
%     figure(FigNums(3));hold on;plot(xinside,yinside,'b.');hold off;
% end

vx=x(edgsInside');
vy=y(edgsInside');

%if FigNums(3)>0
%    figure(FigNums(3));hold on;plot(xinside,yinside,'b.');hold off;
%    hold on;plot(vx,vy,'r-');hold off;
%end

CH=convexHull(DT);
% if FigNums(3)>0
%     figure(FigNums(3));hold on;plot(x(CH), y(CH), 'b');hold off;
% end

edgLengths2 = ( vx(1,:)-vx(2,:) ).^2 + ( vy(1,:)-vy(2,:) ).^2;
meanEdgLength = (mean(edgLengths2.^(-2)))^(-1/4);
disp(['mean particle distance = ' num2str(meanEdgLength)]);
%meanEdgLength = meanEdgLength0;
%disp(['mean edge length =  ' num2str(meanEdgLength)]);
numberOfEntities=length(x);







%% Filter image for second search
% 
% 
% high_limit=cutOff((high_limitFAC/meanEdgLength), 0.03, 1);
% %high_limit=cutOff((high_limitFAC/meanEdgLength)^2, 0.075, 1);
% %high_limit=cutOff(high_limitFAC*numberOfEntities, 0.075, 1);
% 
% h = bp_sample(low_limit2, high_limit, 450);%-0.4*bp_sample(0.0,.01,100);
% %h2=bp_sample(0,.4,25);
% 
% im=filter2fft(h,imSav);
% %im=cutoff(im,0.,2*sqrt(mean(mean(im.^2)) - mean(mean(im))^2));
% % if FigNums(5)>0
% %     figure(FigNums(5));imagesc(im);colorbar
% %     set(FigNums(5),'Position',[610 320 size(imSav0,2)/2. size(imSav0,1)/2.])
% % end
% 
% circleSize100filter = meanEdgLength*.6;
% %thresholdSkew = thresholdSkew2;
% thresholdSkew = cutOff((meanEdgLength-LargeNucDist)/20, 0, 4);
% disp(['threshold Skew = ' num2str(thresholdSkew)]);
% filterMaxNorm
% 
% %% Binarize image and find entities
% figure, imshow(im100)
% threshold    = mean(im100(im100>0)) * thresholdFAC2; %/4;%/10;
% numthreshold = 2;
% preMask   = circle(meanEdgLength*preMaskFAC);
% maskBlackRadius = min(meanEdgLength*maskBlackFAC, maskBlackRadiusMax);
% maskBlack=circle(maskBlackRadius);
% maskWhite = circle(2);
% 
% [x,y,Num,Sum]=entities2(im100,threshold,numthreshold,preMask,maskWhite,maskBlack,[FigNums(7), -1]);


%if FigNums(5)>0
%    figure(FigNums(5));imagesc(cutOff(im,0,1000.0));colorbar;
%    ims=im(YT(channel+1):YB(channel+1),XL(channel+1):XR(channel+1));
%end

%imGrob=filter2(bp_sample(0.0,.01,60),imSav);
%im=im-imGrob;
%if FigNums(5)>0
%    figure(FigNums(5));imagesc(imGrob);colorbar
%    figure(FigNums(5));imagesc(im);colorbar
%end

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



%% Show result and save to file

% if FigNums(1)>0
%     figure(FigNums(1));
%     hold on
%     plot(x,y,'r.');
%     hold off
% else
%     disp('no output png-file written!')
% end

%SAVMAT=[x' y' Num' Sum']; save([outfolder 'outPos' num2str(picnum) '.txt'],'SAVMAT','-ascii', '-tabs');

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

