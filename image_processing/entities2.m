function [X_pivot,Y_pivot,Num,Sum]=entities2(img,treshold,numtreshold,preMask,maskWhite,maskBlack,FigNums)
global imbi imbin image siz num sum x_pivot y_pivot

if nargin<7
    FigNums=[-1 -1];
end
    

X_pivot=NaN;
Y_pivot=NaN;
Num=NaN;
Sum=NaN;
image=img;
siz=size(image);

counter=0;

imbi=image>treshold;

if FigNums(1)>0
    figure(FigNums(1));subplot(2,2,1);hold off;imagesc(imbi);colormap bone
    set(gcf,'Position',[10 10 size(img,2)/2. size(img,1)/2.])
end


%imbi=filter2([1 1 1 1; 1 1 1 1; 1 1 1 1],imbi)~=0;reihentestNR(10,10,2);
if FigNums(2)>0
    figure(FigNums(2));hold off;imagesc(filter2(preMask,imbi));colorbar
end


imbi=filter2(preMask,imbi)>=1;
%reihentestNR(10,10,2);
XxX=2; 
while imbi(XxX,2)~=0
    XxX = XxX + 1;
end
reihentestNR(XxX,2,2);

if FigNums(1)>0
    figure(FigNums(1));subplot(2,2,2);hold off;imagesc(imbi);colormap bone
    %set(gcf,'Position',[300 10 size(img,2)/2. size(img,1)/2.])
end

imbi=~imbi;
if any(any(maskBlack))~=0
    imbi=filter2(maskBlack,imbi)~=0;
end

imbi=~imbi;

if FigNums(1)>0
    figure(FigNums(1));subplot(2,2,3);hold off;imagesc(imbi);colormap bone
    %set(gcf,'Position',[700 10 size(img,2)/2. size(img,1)/2.])
end

if any(any(maskWhite))~=0
    imbi=filter2(maskWhite,imbi)~=0;
    %imbi=~imbi;
end

imbin=int8(imbi);

if FigNums(1)>0
    figure(FigNums(1));subplot(2,2,4);hold off;imagesc(imbin);colormap bone
    %set(gcf,'Position',[700 10 size(img,2)/2. size(img,1)/2.])
end

[yw,xw]=find(imbin==1);
ignored=0;

for i=1:length(yw)
        x=xw(i);
        y=yw(i);
        
        if imbin(y,x)==1
            
            x_pivot=0; y_pivot=0;                       % set pivot zero
    	    num = 0;                                    % number of pixels
	        sum = 0;                                    % gray value sum

                        
            pixelreport(y,x)
            reihentest(y,x)
            
            if num>=numtreshold
                counter=counter+1;
                X_pivot(counter)=x_pivot/sum;
                Y_pivot(counter)=y_pivot/sum;
                Num(counter)=num;
                Sum(counter)=sum;
            else
                ignored=ignored+1;
            end
        end
end
%disp([num2str(ignored) ' elements with less than ' num2str(numtreshold) 'pixels ignored']);
