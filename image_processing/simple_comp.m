close all
for i=[1, 127, 160]
    i;
%data=['/home/lutz/embryo_algo/data/246_interphase_2sec/outPos' num2str(i) '.txt']
% if (i<10)
% movie=['/home/lutz/embryo_algo/movies/246_interphase_2sec/pic-000' num2str(i) '.png']
% end
% if (i>=10 && i< 100)
%     movie=['/home/lutz/embryo_algo/movies/246_interphase_2sec/pic-00' num2str(i) '.png']
% end
% if (i>=100 && i< 1000)
% movie=['/home/lutz/embryo_algo/movies/246_interphase_2sec/pic-0' num2str(i) '.png']
% end
% if (i>=1000 && i< 10000)
% movie=['/home/lutz/embryo_algo/movies/246_interphase_2sec/pic-' num2str(i) '.png']
% end
merge=['/home/lutz/embryo_algo/eval/246_interphase_2sec/pkf' num2str(i) '.txt']
color='green'
switch i
    case 1
        color='red';
    case 127
        color='green';
    case 160
        color='blue';
    otherwise
        color='black';
end
op=load(merge);
figure(2), plot(op(:,1)*0.16,op(:,2),color,'LineWidth',2)
% image=im2double(imread(movie));
% image(1:xx,yy:-1:1)=image(1:xx,1:yy);
% %x=load(data);
% y=load(merge);
% [xx yy]=size(image)
% dummy=image(1:xx,1:yy/2);
% image(1:xx,yy/2:yy)=1;
% %image(find(image==0))=1;
% figure, imshow(image)
hold on

% yc=[]
% [a b]=size(y);
% for j=1:a
%     if(y(j,1) > yy/2)
%         y(j,1)=yy-y(j,1);
%     else
%         y(j,1)=yy-y(j,1);
%     end
% end
% for j=1:a
% % if(y(j,1) > yy/2-60)
% %     yc=[yc; y(j,1) y(j,2)]
% % end
% end
% TRI = delaunay(y(:,1),y(:,2));
% triplot(TRI,y(:,1),y(:,2))
% 
% %plot(x(:,1),x(:,2),'o','MarkerFaceColor','red')
% plot(y(:,1),y(:,2),'o','MarkerFaceColor','green')
% imshow(dummy)
end
legend('t=0','t=252 s','t=318 s')
grid on
xlabel('Distance/\mu m')
ylabel('g(r)')