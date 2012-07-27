function reihentestNR(y,x,n)
global imbi %imbin %image %siz num sum x_pivot y_pivot
%imbi(y,x)
siz=size(imbi);
imbin=int8(imbi);
MAXLISTLEN=1000;
currentListY=zeros(1,MAXLISTLEN);
currentListX=zeros(1,MAXLISTLEN);
currentListNewY=zeros(1,MAXLISTLEN);
currentListNewX=zeros(1,MAXLISTLEN);

currentListLen=1;
currentListY(1)=y;
currentListX(1)=x;


imbin(y,x)=n;

while currentListLen~=0
    currentListNewLen=0;
    for j=1:currentListLen
        yj=currentListY(j);
        xj=currentListX(j);
        
        %search left
        xjnew=xj-1;
        if xj>1 && imbin(yj,xjnew)==0
            imbin(yj,xjnew)=n;
            currentListNewLen=currentListNewLen+1;
            currentListNewY(currentListNewLen)=yj;
            currentListNewX(currentListNewLen)=xjnew;
        end
        
        %search right
        xjnew=xj+1;
        if xj<siz(2) && imbin(yj,xjnew)==0
            imbin(yj,xjnew)=n;
            currentListNewLen=currentListNewLen+1;
            currentListNewY(currentListNewLen)=yj;
            currentListNewX(currentListNewLen)=xjnew;
        end
        
        %search up
        yjnew=yj-1;
        if yj>1 && imbin(yjnew,xj)==0
            imbin(yjnew,xj)=n;
            currentListNewLen=currentListNewLen+1;
            currentListNewY(currentListNewLen)=yjnew;
            currentListNewX(currentListNewLen)=xj;
        end

        %search down
        yjnew=yj+1;
        if yj<siz(1) && imbin(yjnew,xj)==0
            imbin(yjnew,xj)=n;
            currentListNewLen=currentListNewLen+1;
            currentListNewY(currentListNewLen)=yjnew;
            currentListNewX(currentListNewLen)=xj;
        end
    end
    currentListLen=currentListNewLen;
    if currentListNewLen~=0
        currentListY(1:currentListNewLen)=currentListNewY(1:currentListNewLen);
        currentListX(1:currentListNewLen)=currentListNewX(1:currentListNewLen);
    end
end
%figure(206);hold off;imagesc(imbin);colormap bone;colorbar;
imbi=(imbin~=n);
