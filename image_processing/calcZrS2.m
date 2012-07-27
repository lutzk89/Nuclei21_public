function [ Zis, Npartis, Nedgis ] = calcZrS2( x, y, SSym0, r )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
r2=r.*r;
r2max=max(r2);
%pih=pi/2;
n=length(x);
nr=length(r);

NbIndices=cell(n,1);
for ii=1:n
    [NbIndices{ii}, ~, ~]=find(SSym0(:,ii));
end


Zis=zeros(n,nr);
Npartis=zeros(n,nr);
Nedgis=zeros(n,nr);
for ii=1:n
    Nedgi=zeros(1,nr);
    Nparti=zeros(1,nr);
    Zi=zeros(1,nr);
    xi=x(ii);
    yi=y(ii);
    NNbsi=length(NbIndices{ii});
    alphaijs = atan2(y(NbIndices{ii})-yi, x(NbIndices{ii})-xi);
    for kk=1:n
        xk=x(kk);
        yk=y(kk);
        dx=xi-xk;
        dy=yi-yk;
        dr2=dx*dx + dy*dy;
        if dr2 < r2max
            NNbsk = length(NbIndices{kk});
            alphakls = atan2(y(NbIndices{kk})-yk, x(NbIndices{kk})-xk);
            alphaDiffMatrix = alphaijs(ones(NNbsk,1),:) - alphakls(ones(NNbsi,1),:)';
            closeQ = dr2<r2; %Vector, 1 if cond erf.
            Nparti= Nparti+ closeQ;
            Nedgi = Nedgi + closeQ*NNbsi*NNbsk;
            Zi    = Zi+     closeQ*sum(sum(exp(6i*alphaDiffMatrix)));
        end
    end
    Nedgis(ii,:) = Nedgi;
    Npartis(ii,:)= Nparti;
    Zis(ii,:)    = Zi;
end

end
