function [result] = circularSmooth(data,len)
%Written by Gavin Taylor, 2011. MIT License

dataLen = length(data);

%need a catch for short data, but it can work...

newData = [data(end-len+2:end) data data(1:len-1)];
nans = isnan(newData);
newData = smooth(newData', len);
newData(nans) = NaN;
result = newData(len:end-len+1);

