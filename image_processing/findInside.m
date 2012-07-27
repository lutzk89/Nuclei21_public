function [ InsideQ, edgsInside ] = findInside(DT,t)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
n=size(DT.X,1);
edgs=edges(DT);

S=sparse(edgs(:,1),edgs(:,2),1,n,n);
SSym=S+S'+sparse(1:n,1:n,1);

boundaryQ=zeros(n,1);boundaryQ(convexHull(DT))=1;
switch t
    case 0
        boundaryQ2=boundaryQ;
    case 1
        boundaryQ2=SSym*boundaryQ;
    otherwise
        boundaryQ2=SSym^t*boundaryQ;
end
InsideQ=~boundaryQ2;

%xinside=x(InsideQ);
%yinside=y(InsideQ);
%figure(2);hold on;plot(xinside,yinside,'b.');hold off;

boundaryIndices2=find(boundaryQ2~=0);
for i=1:length(boundaryIndices2)
    S(boundaryIndices2(i),:)=0;
    S(:,boundaryIndices2(i))=0;
end
%removeEdgs1=(removeEdgs(:,1)+removeEdgs(:,2)~=0);
%edgsInside=edgs;
[ei,ej,es]=find(S);
edgsInside=[ei ej];
end

