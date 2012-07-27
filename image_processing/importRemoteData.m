function [ out ] = importRemoteData( host, filename, maxDataSize )

if isempty(host)
    out = importdata( filename );
else
    unix(['rsync --quiet -t ' host ':' filename ' /tmp/MatlabTemp.dat']);
    out = importdata('/tmp/MatlabTemp.dat');
    unix('rm /tmp/MatlabTemp.dat');
end

if nargin >= 3 && maxDataSize>0
    DataSize  = size(out,1);
    reduceFac = ceil( DataSize / maxDataSize );
    if reduceFac > 1
        out = out(1:reduceFac:end,:);
    end
end

end

