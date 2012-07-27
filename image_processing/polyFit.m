function ai = polyFit( x, y, degree )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
dataSize = length(x);
deg2p2 = 2*degree+2;
deg1p1 = degree+1;

xPower = ones(deg2p2, dataSize);
for j=2:deg2p2
    xPower(j,:) = xPower(j-1,:).*x;
end

yPower = ones(deg1p1, dataSize);
yPower(1,:) = y;
for j=2:deg1p1
    yPower(j,:) = yPower(j-1,:).*x;
end

sumxPower = sum(xPower,2);
sumyPower = sum(yPower,2);

[i,j] = meshgrid(0:degree,0:degree);
ipj   = i+j;
Mij   = sumxPower(ipj+1);
bj    = sumyPower;
%disp(Mij);
%disp(bj);
ai    = Mij\bj;
ai = ai';
end
