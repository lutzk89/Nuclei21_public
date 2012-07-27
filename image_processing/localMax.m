function [ C, I ] = localMax( array, kernel )
if nargin<2
    convArray = array;
else
    convArray = conv(array,kernel,'same');
end

finiteConvArray = convArray(isfinite(convArray));
meanArray = mean(finiteConvArray);
varArray  = sqrt(var(finiteConvArray))*.01;

Is = zeros(size(array));
for I=(length(array)-2):-1:2
    if convArray(I)>meanArray && convArray(I)>( max( convArray(I-1) , convArray(I+1) )+varArray)
        Is(I)=1;
    end
end

findIs=find(Is);
[~, IsMax] = max(convArray(findIs));

I = findIs(IsMax);
C = array(I);
end

