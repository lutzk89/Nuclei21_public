close all; beep off;
addpath('your folder containing this *.m files')

DataFolder = 'Folder containing the *.tif data files';
% rawfilenames={
% '246_interphase_2sec',  '455_combined',  '456_interphase',  '458_combined',    '620_interphase', '454_combined',	 '456_combined',  '457_combined',    '599_interphase',  '622_time'
% }; % NOTE: no extension (.tif) needed for rawfilenames.

rawfilenames={
'Comma separated list of *.tif files to work on. Have to be in $DataFolder'
    }
for rawfileInd=1:length(rawfilenames)
rawfilename = rawfilenames{rawfileInd};
afilename   = fullfile(DataFolder, [rawfilename '.tif']);
outfolder   = fullfile(DataFolder, rawfilename, filesep);
mkdir(outfolder);

disp(' ');
disp(['filename: ' afilename]);


fileInfo=imfinfo(afilename,'tif');
%include your filename here to cut off timestamps etc before computing
switch rawfilename
    case 'd28'
        [XXX YYY]=ndgrid(1:fileInfo(1).Height,1:fileInfo(1).Width); repaintIndices=find(15*XXX<(YYY-1235)*40); clear XXX YYY;
        repaintValue = 103.;
    case 'd31'
        [XXX YYY]=ndgrid(1:fileInfo(1).Height,1:fileInfo(1).Width); repaintIndices=find(15*XXX<(YYY-1235)*40); clear XXX YYY;
        repaintValue = 110.;
    case 'd38'
        [XXX YYY]=ndgrid(1:fileInfo(1).Height,1:fileInfo(1).Width); repaintIndices=find(XXX-540>(YYY)*60/20); clear XXX YYY;
        repaintValue = 110.;      
    case 'd51'
        [XXX YYY]=ndgrid(1:fileInfo(1).Height,1:fileInfo(1).Width); repaintIndices=find(40-XXX>2*(YYY)); clear XXX YYY;
        repaintValue = 96.;  
    case 'old_data'
        [XXX YYY]=ngrid(1:fileInfo(1).Height,1:fileInfo(1).Width); repaintIndices=ones(YYY,50);
        background=double(imread([DataFolder  'old_bg.tif']));
        background=background/max(max(background));
        sub_bg=1;
        clear XXX YYY;
        repaintValue=26;
    case 'movie1'
        [XXX YYY]=ngrid(1:fileInfo(1).Height,1:fileInfo(1).Width); repaintIndices=ones(YYY,50); clear XXX YYY;
        repaintValue=26;
    case '456_combined' 
        [XXX YYY]=ngrid(1:fileInfo(1).Height,1:fileInfo(1).Width); repaintIndices=ones(YYY,50); clear XXX YYY;
        lutzMask=im2bw(imread([DataFolder '456_mask.png']));
        repaintValue=43;
    case '456_interphase'
        repaintValue=8989;
    case '599_interphase'
        repaintValue=8989;
    otherwise
        repaintValue = -1;
end

FigNums = [1 2 3 100 5 6 100];
%close(FigNums);

colorRange=128;
colors=jet(colorRange);colors(colorRange+1:2000,:)=colors(colorRange*ones(2000-colorRange,1),:);


%% Parameters for the first search 

thresholdSkew1  = 1;
thresholdFAC1   = 1;
low_limit1      = 0.01;

FFTSpectrumTilt = 120.+0*4.;
% If the wrong "image wavelength" is found in figure 2, this parameter tilts (skews) the
% Fourier specturm. 
% Larger values favors larger wavelengths

FFTLocalMaxKernel = 1;5
%FFTLocalMaxKernel = [.25,.5,.25];


%% Parameters for the second search


thresholdFAC2  = .9;
% Threshold for the image binarization
% smaller value: find more nuclei,  connect closeby entities
%  larger value: find fewer nuclei, separate closeby entities

                        

high_limitFAC0 = 4.5;
%high_limitFAC  = 3.;
high_limitFAC  = 4.5;
% larger value: less smoothing, closeby entities are better separated                         
% smaller value: more smoothing, closeby entities are more often 
% considered as one entity


%low_limit2     = .015;
low_limit2     = .035;  
% larger  value: "flatten out" inhomogeneities in illumination
% smaller value: leave image more the way it is


preMaskFAC     = 0.25;   
maskBlackRadiusMax = 18;
maskBlackFAC   = 0.16;   
%  larger value: closeby entities are better separated
% smaller value: closeby entities are more often 
% considered as one entity
                        
                        
thresholdSkew2 = 0.0;   % not very important                        
subtract100fac = 0.6;   % not very important
LargeNucDist   = 20.0;


%% Algorithm Start
picnumMax=length(fileInfo);

%calculate bg
bg=zeros(fileInfo(1).Height,fileInfo(1).Width);
count=0;
for picnum=1:picnumMax
    pic=im2double(imread(afilename,'tif',picnum));
    bg=bg+pic;
    count=count+1;
end
bg=bg/count;
%for picnum=10:10
for picnum=1:picnumMax

disp(' ');
disp(['pic# = ' num2str(picnum)]);



%% Read Image 
imSav0=im2double(imread(afilename,'tif',picnum));
imSav0=imSav0-bg;
%imSav0=im2double(imSav0);

% if repaintValue>=0
%     imSav0(repaintIndices)=repaintValue;
% end

 if repaintValue==26
     imSav0(1:50,:)=0;
 end
 if repaintValue==8989
     imSav0(1:50,1:100)=0;
 end
 [xsize ysize]=size(imSav0);
if repaintValue == 43
%     imSav0(1:xsize,1:100)=mean(mean(imSav0));
%     imSav0(1:50,1:ysize)=mean(mean(imSav0));
%     imSav0(xsize-50:xsize,1:ysize)=mean(mean(imSav0));
%     imSav0(1:xsize,ysize-50:ysize)=mean(mean(imSav0));
    imSav0=imSav0.*lutzMask;
end
% if FigNums(1)>0
%     figure(FigNums(1));hold off;image(imSav0*.29);colormap bone;
%     set(FigNums(1),'Position',[10 630 size(imSav0,2)/2. size(imSav0,1)/2.])
% end

imSav1 = 64*imSav0/max(max(imSav0));
%imSav=imSav-imSavMean;
%imSav=log(imSav+.01);
%imSav     = log(log(imSav0+1)+1);
imSav     = log(imSav0+1);
imSavMean = mean(mean(imSav));%*.65;
imSav     = imSav-imSavMean;
%figure(41);imagesc(imSav);colorbar;

%% Find prevailing wavelength from Fourier Spectrum

meanEdgLength0 = findFrqLog(imSav, 50, FFTSpectrumTilt, FFTLocalMaxKernel, 1., FigNums(2));
%catch error, dunno why this is happening
[corrx corry]=size(meanEdgLength0);
if ~(corry > 0)
    meanEdgLength0=25;
end
if(meanEdgLength0 < 20)
    meanEdgLength0=25;
end
%% Filter image for first search
disp(['prevailing image wave length = ' num2str(meanEdgLength0)]);
high_limit0 = cutOff((high_limitFAC0/meanEdgLength0), 0.032, 1);
h           = bp_sample(low_limit1, high_limit0, 450);
%h           = bp_sample(.05, high_limit0, 450);

im = filter2fft(h,imSav);
if FigNums(3)>0
    figure(FigNums(3));imagesc(im);colorbar;
    %set(FigNums(3),'Position',[610 630 size(imSav0,2)/2. size(imSav0,1)/2.])
end




%% Binarize image and find entities
circleSize100filter = 25;
thresholdSkew = thresholdSkew1;
filterMaxNorm
%threshold    = 0.01;
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
if FigNums(3)>0
    figure(FigNums(3));hold on;plot(xinside,yinside,'b.');hold off;
end

vx=x(edgsInside');
vy=y(edgsInside');

%if FigNums(3)>0
%    figure(FigNums(3));hold on;plot(xinside,yinside,'b.');hold off;
%    hold on;plot(vx,vy,'r-');hold off;
%end

CH=convexHull(DT);
if FigNums(3)>0
    figure(FigNums(3));hold on;plot(x(CH), y(CH), 'b');hold off;
end

edgLengths2 = ( vx(1,:)-vx(2,:) ).^2 + ( vy(1,:)-vy(2,:) ).^2;
meanEdgLength = (mean(edgLengths2.^(-2)))^(-1/4);
disp(['mean particle distance = ' num2str(meanEdgLength)]);
%meanEdgLength = meanEdgLength0;
%disp(['mean edge length =  ' num2str(meanEdgLength)]);
numberOfEntities=length(x);






if FigNums(6)>0
    figure(FigNums(6));subplot(2,2,3);imagesc(im100);colorbar
end

%% Filter image for second search


high_limit=cutOff((high_limitFAC/meanEdgLength), 0.03, 1);
%high_limit=cutOff((high_limitFAC/meanEdgLength)^2, 0.075, 1);
%high_limit=cutOff(high_limitFAC*numberOfEntities, 0.075, 1);

h = bp_sample(low_limit2, high_limit, 450);%-0.4*bp_sample(0.0,.01,100);
%h2=bp_sample(0,.4,25);

im=filter2fft(h,imSav);
%im=cutoff(im,0.,2*sqrt(mean(mean(im.^2)) - mean(mean(im))^2));
if FigNums(5)>0
    figure(FigNums(5));imagesc(im);colorbar
    %set(FigNums(5),'Position',[610 320 size(imSav0,2)/2. size(imSav0,1)/2.])
end

circleSize100filter = meanEdgLength*.6;
%thresholdSkew = thresholdSkew2;
thresholdSkew = cutOff((meanEdgLength-LargeNucDist)/20, 0, 4);
disp(['threshold Skew = ' num2str(thresholdSkew)]);
filterMaxNorm

%% Binarize image and find entities
threshold    = mean(im100(im100>0)) * thresholdFAC2; %/4;%/10;
numthreshold = 2;
preMask   = circle(meanEdgLength*preMaskFAC);
maskBlackRadius = min(meanEdgLength*maskBlackFAC, maskBlackRadiusMax);
maskBlack=circle(maskBlackRadius);
maskWhite = circle(2);

[x,y,Num,Sum]=entities2(im100,threshold,numthreshold,preMask,maskWhite,maskBlack,[FigNums(7), -1]);


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

if FigNums(1)>0
    figure(FigNums(1));
    hold on
    plot(x,y,'r.');
    hold off
    set(gcf, 'PaperPositionMode', 'auto');
    print(gcf,'-dpng',[outfolder 'pic' num2str(picnum) '.png'],'-r200')
else
    disp('no output png-file written!')
end

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

end
if repaintValue==26
    background=background/double(number);
    imwrite(background,[DataFolder + '/background_old_data.tif']);
end
