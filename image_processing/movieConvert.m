obj = VideoReader('/data/hades/data/ulrich/Coworkers/Vincenzo/slinky4512K_Stream.avi')

for i=1:2205
    frame = read(obj, i); 
    imwrite(frame, ['/data/hades/data/ulrich/Coworkers/Vincenzo/slinky/slinky' num2str4(i) '.png'],'png');
end
%myVideo = VideoWriter('/data/hades/data/ulrich/Coworkers/Vincenzo/myfile.avi');
%open(myVideo);
%for i=1:62; frame = read(obj, i); imwrite(frame, ['~/Downloads/mov1/mov1_' num2str(i) '.png'],'png');end