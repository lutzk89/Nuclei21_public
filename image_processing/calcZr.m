function [ Zis, Npartis, Nedgis ] = calcZr( DT, r )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
r2=r*r;
pih=pi/2;
x=DT.X(:,1)';
y=DT.X(:,2)';
n=length(x);

edgs=edges(DT);
S=sparse(edgs(:,1),edgs(:,2),1,n,n);
SSym0=S+S';

NbIndices=cell(n,1);
for ii=1:n
    [NbIndices{ii}, ~, ~]=find(SSym0(:,ii));
end


Zis=zeros(n,1);
Npartis=zeros(n,1);
Nedgis=zeros(n,1);
for ii=1:n
    Nedgi=0;
    Nparti=0;
    Zi=0;
    xi=x(ii);
    yi=y(ii);
    NNbsi=length(NbIndices{ii});
    alphaijs=atan2(y(NbIndices{ii})-yi,x(NbIndices{ii})-xi);
    for kk=1:n
        xk=x(kk);
        yk=y(kk);
        dx=xi-xk;
        dy=yi-yk;
        if dx*dx + dy*dy < r2
            NNbsk=length(NbIndices{kk});
            alphakls=atan2(y(NbIndices{kk})-yk,x(NbIndices{kk})-xk);
            Nparti= Nparti+ 1;
            Nedgi = Nedgi + NNbsi*NNbsk;
            alphaijs/pi*180
            alphakls/pi*180
            kk
            NbIndices{kk}
            alphaDiffMatrix=alphaijs(ones(NNbsk,1),:)-alphakls(ones(NNbsi,1),:)';
            %alphaDiffMatrix=mod(alphaDiffMatrix+pih,pi)-pih;
            alphaDiffMatrix/pi*180;
            exp(6i*alphaDiffMatrix)
            Zi=Zi+sum(sum(exp(6i*alphaDiffMatrix)));
        end
    end
    Nedgis(ii) =Nedgi;
    Npartis(ii)=Nparti;
    Zis(ii)=Zi;
end

end
