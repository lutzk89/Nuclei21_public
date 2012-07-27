function out=circle2x(s1,w1,s2,w2)
rsmax=round(max(s1,s2));
rsmaxh=rsmax/2.;
out=zeros(rsmax);
ss1=s1*s1/4;
ss2=s2*s2/4;
for i=1:rsmax
    for j=1:rsmax
        if (i-.5-rsmaxh)^2+(j-.5-rsmaxh)^2<ss1
            out(i,j)=out(i,j)+w1;
        end
        if (i-.5-rsmaxh)^2+(j-.5-rsmaxh)^2<ss2
            out(i,j)=out(i,j)+w2;
        end
    end
end
