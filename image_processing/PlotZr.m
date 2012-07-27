edgs=edges(DT);S=sparse(edgs(:,1),edgs(:,2),1,n,n);SSym0=S+S';
InsideIndices=find(InsideQ);
SSym0Inside=SSym0(InsideIndices,InsideIndices);

[Zr,Np,Ne]=calcZrS2( xInside, yInside, SSym0Inside , r_Zr);
ZrNpNe{picnum,1}=Zr;
ZrNpNe{picnum,2}=Np;
ZrNpNe{picnum,3}=Ne;
%ZrAbs(picnum,:)=abs(sum(Zr./Ne,1))/length(xInside);
ZrAbs(picnum,:)=abs(sum(Zr,1)./sum(Ne,1));

subplot(1,2,1);image(imSav1);colormap bone;
axis([0 size(imSav1,2) 0 size(imSav1,1)]);axis equal;axis image;
hold on;
plot(xInside,yInside,'r.');%,'MarkerSize',3);
plot(x(CH), y(CH), 'b', 'LineWidth',1);
hold off;
subplot(1,2,2);
loglog(r_Zr,ZrAbs(picnum,:),'.-');
xlabel('r');ylabel('Z(r)');title(['t=' num2str4(picnum)]);
axis square;axis([20 604 0.0005 .5]);%axis([20 1000 0.0001 .5])
