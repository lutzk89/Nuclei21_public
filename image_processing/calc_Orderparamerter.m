n  = size(pos,1);
nP = size(posP,1);
InsideQ = zeros(size(posP,1),1);
InsideQ(1:n) = 1;
InsideIndices = 1:n;

x = posP(:,1)';
y = posP(:,2)';

DT = DelaunayTri(x',y');
[V,R] = voronoiDiagram(DT);
%[ InsideQ, edgsInside ] = findInsideCH(DT,2);
%InsideIndices=find(InsideQ);

edgs=edges(DT);
S=sparse(edgs(:,1),edgs(:,2),1,nP,nP);
SSym0=S+S';

[Zr, Np, Ne] = calcZrS2( x, y, SSym0, r_Zr);

Zr = Zr(InsideIndices,:);
Np = Np(InsideIndices,:);
Ne = Ne(InsideIndices,:);


NNeighbors=full(sum(SSym0,1))';
Inside6=(InsideQ & NNeighbors==6);
Inside5=(InsideQ & NNeighbors<=5);
Inside7=(InsideQ & NNeighbors>=7);


absZr=abs(Zr./Ne);
%MeanAbsZr=mean(absZr);
MeanAbsZr=abs(sum(Zr,1,'double')./sum(Ne,1,'double'));

