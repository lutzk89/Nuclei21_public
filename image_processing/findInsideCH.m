function [ InsideQ, edgsInside ] = findInsideCH(DT,t)
t2=t*t;
X=DT.X; n=size(X,1);
edgs=edges(DT);
CH=convexHull(DT);

InsideQ=ones(n,1);

for i=1:length(CH)-1
    b=X(CH(i),:);
    c=X(CH(i+1),:);
    cb=c-b;
    Xb=X-b(ones(n,1),:);

    cb2=sum(cb.^2);
    Xb2=sum(Xb.^2,2);
    Xbcb=sum(Xb.*cb(ones(n,1),:),2);

    d2=Xb2-(Xbcb.^2)/cb2;
    InsideQ=InsideQ & d2>t2;
end


boundaryQ2=~InsideQ;
boundaryIndices2=find(boundaryQ2~=0);

S=sparse(edgs(:,1),edgs(:,2),1,n,n);
for i=1:length(boundaryIndices2)
    S(boundaryIndices2(i),:)=0;
    S(:,boundaryIndices2(i))=0;
end
[ei,ej,es]=find(S);
edgsInside=[ei ej];
end

