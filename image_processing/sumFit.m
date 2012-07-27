function [ai, yFit] = sumFit( x, y, powerArray )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
degree = length(powerArray);

sumy = zeros(degree, 1);
for j=1:degree
    sumy(j,1) = sum(y.*(x.^powerArray(j)));
end

Mij = zeros(degree, degree);
for i=1:degree
    for j=1:degree
        Mij(i,j) = sum((x.^powerArray(i)) .* (x.^powerArray(j)));
    end
end

ai = Mij\sumy;
ai = ai';

if nargout >= 2
    yFit = zeros(size(x));
    for j=1:degree
        yFit = yFit + ai(j)*x.^powerArray(j);
    end
end

end
