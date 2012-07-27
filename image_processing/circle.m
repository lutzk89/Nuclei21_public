function out=circle(s)
rs=round(s);
out=zeros(rs);
sh=rs/2;
ss=s*s/4;
for i=1:rs
    for j=1:rs
        if (i-.5-sh)^2+(j-.5-sh)^2<ss
            out(i,j)=1;
        end
    end
end
