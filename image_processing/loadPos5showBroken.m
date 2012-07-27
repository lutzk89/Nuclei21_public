clear
hostname  = ''; cutInds = []'; offset = 0;
folder0='/home/lutz/Nuclei21dir/';

folderAdd    = 'run/';
outfolderAdd = 'movLR/';

colorRange    = 128;
colorRangeEXT = 2000;
colors=jet(colorRange);colors(colorRange+1:colorRangeEXT,:)=colors(colorRange*ones(colorRangeEXT-colorRange,1),:);
colors(colorRangeEXT+1,:) = [1 1 1];

%relevant
showConn_switch  = 1; ConnColor  = [1 0 0];
showConn0_switch = 0; Conn0Color = [1 0 0];
showOrder_switch = 1;
showBoundary_switch = 1;
orderParameter_switch = 1;
write_switch = 1;
simplePlot_overwriteSwitch = 0;
%relevant
imageSizeX = 800;
imageBitmapSize = 300;
markerSize = 3;
%markerSize = .2;
%markerSize = .5;
lineWidth  = 2.2;
outResolution = 150;
countMax = 8000;
shearRateFac = 0.10;


vInd      = 1;
stress_ij = [2 1];
strain_ij = [1 2];

if isempty(hostname)
    remoteFolder0 = folder0;
else
    remoteFolder0 = fullfile('/clusterdata/ulrich/',folder0(25:end));
end
outfolder = fullfile(folder0, outfolderAdd);
fid = fopen(fullfile(folder0, 'Settings.txt')); C = textscan(fid, '%s%s%s%s'); fclose(fid);
Settings = cell2struct(C{2},C{1});
dtPos   = str2num(Settings.dt_pos);
dtForce = max(0.001,str2num(Settings.dt_force));

SettingsBS    = cell2struct([C{2} C{3} C{4}],C{1});
[BS1 BS2 BS3] = SettingsBS.BoxSize; 
if isempty(BS3);
    dimension = 2;
else
    dimension = 3; 
end


stressInd = (stress_ij(1)-1) * dimension + stress_ij(2);
strainInd = (strain_ij(1)-1) * dimension + strain_ij(2);

BOUNDARYpos   = importRemoteData(hostname, fullfile(remoteFolder0, folderAdd, 'BOUNDARYpos.txt'));
BOUNDARYforce = importRemoteData(hostname, fullfile(remoteFolder0, folderAdd, 'BOUNDARYf.txt'));
if exist([folder0 '../BOUNDARYm' '.txt'],'file')
    BOUNDARYshear = importdata([folder0 '../BOUNDARYm' '.txt']);
elseif exist([folder0 '../BOUNDARYo' '.txt'],'file') 
    BOUNDARYshear = importdata([folder0 '../BOUNDARYo' '.txt']);
else
    BOUNDARYshear = zeros(1,3*dimension^2 + 1);
end

shearRate = BOUNDARYshear(1,strainInd);
countMax  = min(countMax,round(size(BOUNDARYpos,1)/dtPos*dtForce));
disp([' -- countMax: ' num2str(countMax) ' -- ']);
L0      = abs(BOUNDARYpos(1,[1 dimension+2])');
L       = L0;
Lh      = L/2;
invRim  = [1; 1];
rim     = L .* invRim;
rimPlot = [3*L(1) 3*L(1) L(2) L(2)]/12;
%rimPlot = [L(1) L(1) L(2) L(2)]/50;
times   = (0:countMax)*dtPos;
r_Zr    = [0.0001 0.5:0.5:0.5*min(rim)];
MeanAbsZrs  = zeros(countMax+1,length(r_Zr));
VelVariance = [];

figure(1); clf;
set(1,'Position',[10 10 imageSizeX round(imageSizeX * (L0(2)+rimPlot(4)-rimPlot(3))/(L0(1)+rimPlot(2)-rimPlot(1)))]);

if showConn0_switch > 0
    conn0 = importRemoteData(hostname, fullfile(remoteFolder0, folderAdd, ['OUTc' num2str(0) '.txt']));
    conn0Len = size(conn0,1);
else
    conn0Len = 0;
end

for count = 0:1:countMax-1
    %fCount = round(count*dtPos/dtForce)+1;
    fCount=1;
    disp(['count: ' num2str(count) '     fCount: ' num2str(fCount)]);
    LEshift    = reshape(BOUNDARYpos(fCount,1:dimension^2),[dimension dimension])';
    invLEshift = inv(LEshift);
    L       = abs(BOUNDARYpos(fCount,[1 4])');
    Lh      = L/2;

    pos0   = importRemoteData(hostname, fullfile(remoteFolder0, folderAdd, ['OUTx' num2str(count) '.txt']));
    vel0   = importRemoteData(hostname, fullfile(remoteFolder0, folderAdd, ['OUTv' num2str(count) '.txt']));
    stressFile   =                      fullfile(remoteFolder0, folderAdd, ['OUTs' num2str(count) '.txt']);
    velFluctFile =                      fullfile(remoteFolder0, folderAdd, ['OUTvf' num2str(count) '.txt']);
    %if exist(stressFile,'file')
        stress = importRemoteData(hostname, stressFile);
    %else
    %    stress = zeros(size(pos0));
    %end
    %if exist(stressFile,'file') vel0Fluct = importRemoteData(hostname, velFluctFile);
    %else                        vel0Fluct = zeros(size(vel0)); end
    
    vel0Fluct = zeros(size(vel0));

    
    pos0(:,2) = pos0(:,2) + offset*Lh(2);
    VelVariance(count+1,:) = [var(vel0(:,1)) var(vel0(:,2))];
    
    
    
    set(1,'Name',['time ' num2str(times(count+1))],'NumberTitle','off');
        
    if simplePlot_overwriteSwitch > 0
        maxv = +shearRateFac*shearRate;
        minv = -shearRateFac*shearRate;
        maxv =  0.000001+2*sqrt(VelVariance(count+1,vInd));
        minv = -0.000001-2*sqrt(VelVariance(count+1,vInd));
        %plot(pos0(1,:),pos0(2,:),'ko','MarkerFaceColor','k','MarkerSize',markerSize);
        pos0(cutInds,:) = [];
        
        
        switch simplePlot_overwriteSwitch
            case 1
                scatter(pos0(:,1), pos0(:,2), 15*markerSize, vel0(:,vInd),'filled');
                axis([-rimPlot(1) L(1)+rimPlot(2) -rimPlot(3) L(2)+rimPlot(4)]);
                set(gca,'CLimMode','manual','CLim',[minv maxv]);
                
            case 2
                pos = (invLEshift*pos0')';
                scatter(pos(:,1), pos(:,2), 15*markerSize, vel0(:,vInd), 'filled');
                axis([-.1 1.1 -.1 1.1]);
                set(gca,'CLimMode','manual','CLim',[minv maxv]);
                
            case 3
                figure(11);clf;
                scatter(pos0(:,2), vel0(:,1),'.');
                if write_switch<0
                    print(gcf,'-dpng',[outfolder 'profile' num2str(count) '.png'],['-r' num2str(outResolution)]);
                end
                
            case 4
                maxStress = .005000;
                minStress = -maxStress;
                %zdata0 = -stress(:,stressInd) - stress(:,stressInd+4) + vel0Fluct(:,stress_ij(1)).*vel0Fluct(:,stress_ij(2));
                zdata0 = vel0Fluct(:,stress_ij(1)).*vel0Fluct(:,stress_ij(2));
                
                scatter(pos0(:,1), pos0(:,2), 15*markerSize, zdata0 , 'filled');
                axis([-rimPlot(1) L(1)+rimPlot(2) -rimPlot(3) L(2)+rimPlot(4)]);
                set(gca,'CLimMode','manual','CLim',[minStress maxStress]);
                
            case 1005
                maxStress = .3;
                minStress = -.1;%maxStress;
                zdata0 = - ( stress(:,stressInd) + stress(:,stressInd+4) );
                zdata = colorRange * ( zdata0 - minStress ) / ( maxStress - minStress );
                
                imageBitmapSizeX = round( imageBitmapSize * (1 + (rimPlot(1)+rimPlot(2))/L0(1) ) );
                imageBitmapSizeY = round( imageBitmapSize * (1 + (rimPlot(3)+rimPlot(4))/L0(2) ) );
                
                imageIndices = [round( imageBitmapSize * (pos0(:,1)+rimPlot(1)) / L0(1))  ,  round( imageBitmapSize * (pos0(:,2)+rimPlot(3)) / L0(2))];
                imageBitmap  = accumarray(imageIndices , zdata , [imageBitmapSizeX imageBitmapSizeY], @mean, colorRangeEXT+1);
                
                image(imageBitmap');set(gca,'YDir','normal');colormap(colors);
                %colorbar;
                %scatter(pos0(:,1), pos0(:,2), 15*markerSize, -stress(:,stressInd) - stress(:,stressInd+4), 'filled');
                %axis([-rimPlot(1) L(1)+rimPlot(2) -rimPlot(3) L(2)+rimPlot(4)]);
                %set(gca,'CLimMode','manual','CLim',[minStress maxStress]);

            case 1006
                maxz = 10;
                minz = -1;%maxStress;
                zdata0 = pos0(:,end);
                zdata = colorRange * ( zdata0 - minz ) / ( maxz - minz );
                
                imageBitmapSizeX = round( imageBitmapSize * (1 + (rimPlot(1)+rimPlot(2))/L0(1) ) );
                imageBitmapSizeY = round( imageBitmapSize * (1 + (rimPlot(3)+rimPlot(4))/L0(2) ) );
                
                imageIndices = [round( imageBitmapSize * (pos0(:,1)+rimPlot(1)) / L0(1))  ,  round( imageBitmapSize * (pos0(:,2)+rimPlot(3)) / L0(2))];
                imageBitmap  = accumarray(imageIndices , zdata , [imageBitmapSizeX imageBitmapSizeY], @mean, colorRangeEXT+1);
                
                image(imageBitmap');set(gca,'YDir','normal');colormap(colors);
                %colorbar;
                %scatter(pos0(:,1), pos0(:,2), 15*markerSize, -stress(:,stressInd) - stress(:,stressInd+4), 'filled');
                %axis([-rimPlot(1) L(1)+rimPlot(2) -rimPlot(3) L(2)+rimPlot(4)]);
                
            case 2006
                pos0(:,2) = L0(2) - pos0(:,2);
                maxz = 10; minz = -1;
                %maxz = 6; minz = -6;
                zdata0 = pos0(:,end);
                
                maxStress = .3;
                minStress = -.1;
                zdata1  = - ( stress(:,stressInd) + stress(:,stressInd+4) );
                zdata1s = colorRange * ( zdata1 - minStress ) / ( maxStress - minStress );

                
                imageBitmapSizeX = round( imageBitmapSize * (1 + 0*(rimPlot(1)+rimPlot(2))/L0(1) ) );
                imageBitmapSizeY = round( imageBitmapSize * (1 + 0*(rimPlot(3)+rimPlot(4))/L0(2) ) );
                
                imageIndices = ceil([ imageBitmapSize * (pos0(:,1)+0*rimPlot(1)) / L0(1)  , imageBitmapSize * (pos0(:,2)+0*rimPlot(3)) / L0(2)] );
                goodIndices  = all( 0<imageIndices, 2) & imageIndices(:,1)<=imageBitmapSizeX & imageIndices(:,2)<=imageBitmapSizeY;
                imageBitmap0  = accumarray(imageIndices(goodIndices,:) , zdata0(goodIndices,:)  , [imageBitmapSizeX imageBitmapSizeY], @mean, NaN);
                imageBitmap1  = accumarray(imageIndices(goodIndices,:) , zdata1s(goodIndices,:) , [imageBitmapSizeX imageBitmapSizeY], @mean, NaN);
                
                %goodIndices = ~isnan(imageBitmap0);
                
                
                [meshX meshY] = meshgrid( (1:imageBitmapSizeX)/imageBitmapSizeX*L0(1), (1:imageBitmapSizeY)/imageBitmapSizeY*L0(2) );
                %surf(meshX(goodIndices), meshY(goodIndices), imageBitmap0(goodIndices), 'LineStyle','none');
                surf(meshX, meshY, imageBitmap0, 'LineStyle','none');
                axis([0 L0(1) 0 L0(2) minz maxz minz maxz])

                %surf(imageBitmap0,imageBitmap1,'LineStyle','none');

            otherwise
                disp('get lost... sucker');
                
        end
        
        
    else
        
        %Larray = L*ones(1,size(pos0,2));
    
    
        invPos0       = (invLEshift*pos0')';
        PeriodicCount = floor( invPos0 );
        
        invPos  = invPos0 - PeriodicCount;
        pos     = (LEshift*invPos')';
        
        invPosP = invPosPeriodicLE(invPos, [0.5 0.5]);
        posP    = (LEshift*invPosP')';
        
        %scatter(posP(:,1),posP(:,2),'.');
        %i feel, that the middle of the map should not be empty. but possibly i
        %am causing bugs doing this
        posP=[pos0; posP];

        if showConn_switch > 0
            conn = importRemoteData(hostname, fullfile(remoteFolder0, folderAdd, ['OUTc' num2str(count) '.txt']));
            connLen = size(conn,1);
        else
            connLen = 0;
        end


        if orderParameter_switch > 0
            calc_Orderparamerter;
            MeanAbsZrs(count+1,:) = MeanAbsZr;
        end





        if showBoundary_switch > 0
            LEboundary = LEshift;
            boundary = (LEboundary * [0 1 1 0 0;
                                      0 0 1 1 0])';
%                         LEboundary*[1 0];
%                         LEboundary*[1 1];
%                         LEboundary*[0 1];
%                         LEboundary*[0 0]];
%             boundary = [[0 0]*LEboundary;
%                         [1 0]*LEboundary;
%                         [1 1]*LEboundary;
%                         [0 1]*LEboundary;
%                         [0 0]*LEboundary;];

            plot(boundary(:,1),boundary(:,2),'b');
            hold on;
        end
        
        
        if showConn0_switch > 0 && conn0Len > 0
            connBSparse = sparse(conn0(:,1),conn0(:,2),1,size(pos0,1),size(pos0,1)) - sparse(conn(:,1),conn(:,2),1,size(pos0,1),size(pos0,1));
            [connB1 connB2] = find(connBSparse);
            connBLen = length(connB2);
            Lhs   = ones(connBLen,1)*Lh';
            Ls    = ones(connBLen,1)*L';
            P1    = pos(connB1,:);
            P2    = pos(connB2,:);
            dP    = P2 - P1;
            Periodic   = round(dP ./ Ls);
            NPeriodicQ = all(Periodic==0,2);


            plot([P1(NPeriodicQ,1) P2(NPeriodicQ,1)]',[P1(NPeriodicQ,2) P2(NPeriodicQ,2)]','Color',Conn0Color,'LineWidth', showConn0_switch*lineWidth);
            hold on;

            dPs = dP - (LEshift*Periodic')';
            for ci = (find(~NPeriodicQ)')
                P1a = P1(ci,:);
                P2a = P1a + dPs(ci,:);

                P2b = P2(ci,:);
                P1b = P2b - dPs(ci,:);

                plot([P1a(1,1) P2a(1,1); P1b(1,1) P2b(1,1)]',[P1a(1,2) P2a(1,2); P1b(1,2) P2b(1,2)]','Color',0.4*Conn0Color + 0.6*[1 1 1],'LineWidth',lineWidth);
            end


        end

        if showConn_switch > 0 && connLen > 0
            Lhs   = ones(connLen,1)*Lh';
            Ls    = ones(connLen,1)*L';
            P1    = pos(conn(:,1),:);
            P2    = pos(conn(:,2),:);
            dP    = P2 - P1;
            Periodic   = round(dP ./ Ls);
            NPeriodicQ = all(Periodic==0,2);


            plot([P1(NPeriodicQ,1) P2(NPeriodicQ,1)]',[P1(NPeriodicQ,2) P2(NPeriodicQ,2)]','Color',ConnColor,'LineWidth',lineWidth);
            hold on;

            dPs = dP - (LEshift*Periodic')';
            for ci = (find(~NPeriodicQ)')
                P1a = P1(ci,:);
                P2a = P1a + dPs(ci,:);

                P2b = P2(ci,:);
                P1b = P2b - dPs(ci,:);

                plot([P1a(1,1) P2a(1,1); P1b(1,1) P2b(1,1)]',[P1a(1,2) P2a(1,2); P1b(1,2) P2b(1,2)]','Color',0.4*ConnColor + 0.6*[1 1 1],'LineWidth',lineWidth);
            end


        end
        hold on;

        plot( pos(:,1), pos(:,2),'ko','MarkerFaceColor','k','MarkerSize',markerSize);
        plot(posP(:,1),posP(:,2),'o', 'Color',[.5 .5 .5],'MarkerSize',markerSize);

        if showOrder_switch>0
            for Ind=1:n
                patch(V(R{Ind},1),V(R{Ind},2),[1 0 0],'FaceColor',colors(floor(colorRange*absZr(Ind))+1,:));
            end
        end

        hold off;
        axis equal;
        
        axis([-rimPlot(1) L0(1)+rimPlot(2) -rimPlot(3) L0(2)+rimPlot(4)]);
        %axis([-4 8 -4 8]);
        




        if showOrder_switch>0
            figure, loglog(r_Zr,MeanAbsZr)%;axis([0.4 max(r_Zr)+.1 0.02 1])
        end
    end
    
    if gcf ~= 1; figure(1); end
    hold off;
	if simplePlot_overwriteSwitch < 2000
        axis equal;
    end
    
    if simplePlot_overwriteSwitch < 1000
        axis([-rimPlot(1) L0(1)+rimPlot(2) -rimPlot(3) L0(2)+rimPlot(4)]);
    end
    
    if write_switch>0
        if showConn_switch > 0 && simplePlot_overwriteSwitch == 0
            print(gcf,'-dpng',[outfolder 'pcon' num2str(count) '.png'],['-r' num2str(outResolution)]);
        else
            print(gcf,'-dpng',[outfolder 'particles' num2str(count) '.png'],['-r' num2str(outResolution)]);
        end
    end
    
end


if orderParameter_switch > 0
    %figure(10);plot(MeanAbsZrs); axis([-1 countMax+1 -.01 1.01])
    figure(11);semilogy(times,MeanAbsZrs(:,2:size(MeanAbsZrs,2)));axis([-1 times(countMax)+1 0.005 1])
end


if write_switch>1
    print(gcf,'-dpng', fullfile(folder0, outfolderAdd, 'Zrs.eps'       ));
    save(              fullfile(folder0, outfolderAdd, 'MeanAbsZrs.mat'),'r_Zr','MeanAbsZrs');
end

%figure(2); hold off; plot(VelVariance(:,1)); hold on; plot(100*VelVariance(:,2),'r');
forceIndex = 4;
%figure(3); hold off; plot(BOUNDARYpos(:,2),BOUNDARYforce(:,forceIndex)+BOUNDARYforce(:,forceIndex+4)); hold on; %plot(BOUNDARYpos(:,2),100*BOUNDARYforce(:,5),'r');
