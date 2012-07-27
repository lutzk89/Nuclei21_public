function [ out ] = importClusterData( host, filename )

if host==''
    out = importdata( filename );
else
    unix(['rsync -t ' host ':' filename ' /tmp/ulrich/MatlabTemp.dat']);
    out = importdata('/tmp/ulrich/MatlabTemp.dat');
    unix('rm /tmp/ulrich/MatlabTemp.dat');
end

end

