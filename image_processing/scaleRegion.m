function [ outArray ] = scaleRegion( inArray, factor )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
inArrayMean     = mean(inArray);
inArrayCentered = inArray - inArrayMean;
outArray        = inArrayCentered * factor + inArrayMean;

end

