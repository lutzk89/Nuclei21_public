% The following variables have to be defined:
% x,y:     particle positions
% DT:      Delaunay Triangulation
% InsideQ: boolean vector, stating which partices are considered for 
%          data analysis
% meanEdgLength: typical length between particles
% colors:  color scheme for the voronoi patches: colorRange x 3 matrix

[V,R] = voronoiDiagram(DT);
n=length(x);
edgs=edges(DT);
S=sparse(edgs(:,1),edgs(:,2),1,n,n);
SSym0=S+S';

xInside=x(InsideQ);
yInside=y(InsideQ);
InsideIndices=find(InsideQ);
SSym0Inside=SSym0(InsideIndices,InsideIndices);
%markerSize=meanEdgLength/2.
markerSize=10;


%%%%%%%%%%%%% Plot 1 %%%%%%%%%%%%%%%
subplot(2,2,1);image(imSav1);colormap bone;
axis([0 size(imSav1,2) 0 size(imSav1,1)]);axis equal;axis image;
hold on;
plot(xInside,yInside,'r.');%,'MarkerSize',3);
plot(x(CH), y(CH), 'b', 'LineWidth',1);
hold off;



%%%%%%%%%%%%% Plot 2 %%%%%%%%%%%%%%%
subplot(2,2,2);image(imSav1);colormap bone;

densities=zeros(n,1);
for ii=1:length(InsideIndices)
    Ind=InsideIndices(ii);
    densities(Ind)=1/polyarea(V(R{Ind},1),V(R{Ind},2));
end

MeanDensity=8e-04;
densities=cutOff(densities/MeanDensity,0,2);

axis([0 size(imSav1,2) 0 size(imSav1,1)]);axis equal;axis image;
hold on;
for ii=1:length(InsideIndices)
    Ind=InsideIndices(ii);
    patch(V(R{Ind},1),V(R{Ind},2),[1 0 0],'FaceColor',colors(round(colorRange*.5*densities(Ind)),:));
end
hold off;



%%%%%%%%%%%%% Plot 3 %%%%%%%%%%%%%%%
subplot(2,2,3);
axis([0 size(imSav1,2) 0 size(imSav1,1)]);axis equal;axis image;axis ij;
triplot(DT,'Color','k');

NNeighbors=full(sum(SSym0,1))';
Inside6=(InsideQ & NNeighbors==6);
Inside5=(InsideQ & NNeighbors<=5);
Inside7=(InsideQ & NNeighbors>=7);

hold on;
plot(x(Inside6),y(Inside6),'.','MarkerSize',markerSize,'Color',[1 0 0]);
plot(x(Inside5),y(Inside5),'.','MarkerSize',markerSize,'Color',[0 .7 .3]);
plot(x(Inside7),y(Inside7),'.','MarkerSize',markerSize,'Color',[0 .3 .7]);
hold off;
%axis equal;axis ij;
axis([0 size(imSav1,2) 0 size(imSav1,1)]);axis image;axis equal;axis ij;




%%%%%%%%%%%%% Plot 4 %%%%%%%%%%%%%%%
subplot(2,2,4);
%axis([0 size(imSav1,2) 0 size(imSav1,1)]);axis equal;axis image;axis ij;
plot(x(CH), y(CH), 'k', 'LineWidth',1);



[Zr,Np,Ne]=calcZrS( xInside, yInside, SSym0Inside , 1.5*meanEdgLength);
absZr=abs(Zr./Ne);
MeanAbsZr=mean(absZr);
absZr=cutOff(2*absZr,0,1);
%absZr=cutOff(log10(absZr)+5,0,5)/5;
%absZr=cutOff(absZr/MeanAbsZr,0,2);
hold on;
for ii=1:length(InsideIndices)
    Ind=InsideIndices(ii);
    patch(V(R{Ind},1),V(R{Ind},2),[1 0 0],'FaceColor',colors(1+round(128*absZr(ii)),:));
end
hold off;
axis([0 size(imSav1,2) 0 size(imSav1,1)]);axis image;axis equal;axis ij;

% hold on;
% for ii=1:length(InsideIndices)
%     Ind=InsideIndices(ii);
%     patch(V(R{Ind},1),V(R{Ind},2),[1 0 0],'FaceColor',colors(round(64*areas(Ind)),:));
% end
% hold off;
% axis([0 size(imSav1,2) 0 size(imSav1,1)]);axis equal;axis image;