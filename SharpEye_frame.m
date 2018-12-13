% read in Kelvin Hughes Sharp Eye radar data, generate Frequency 1 and 
% Frequency 2 scans and save to disk, ECEF (Earth centric Earth fixed) 
% registered with ship data (lat, lon, heading)

clear all
close all

% get the parameters
load('parameters.mat');

% open a file               
% d = dir([datapath '*.dat']);
datafile = uigetfile([datapath '*.dat']);

% load the file index
% di = dir([datapath '*.mat']);
load([datapath datafile(1:6) '_index.mat']);

%load in the ShipData file for the date of interest
[sdnum, slat, slon, scog, ssog, shead] = load_shipData(datapath(18:25));

for ns=1:length(index)
    tic
    [scandnum,...
     f1time,f1mag,f1phi,...
     f2time,f2mag,f2phi,...
     fileheader,...
     minf1turning,maxf1turning,...
     minf2turning,maxf2turning,...
     squint,...
     r,theta,...
     lonp,latp] = read_SharpEye_sector([datapath datafile], index(ns,2));
    toc
    
    % find range and azimuth indexes to process, defines a wedge, from parameters
    r1 = find(r>=rmin,1,'first');
    r2 = find(r>=rmax,1,'first');
    [val, az1] = min(abs(theta-azmin*pi/180));
    [val, az2] = min(abs(theta-azmax*pi/180));
    
    % if we're straddling the scan boundary, reduce the wedge size
    if az2<az1
       taz = 4096-az1;
       naz = floor(taz./squint);
       az2 = az1 + (naz.*squint);
    end
    
    theta = unwrap(theta);
    mtheta = pi/2 - theta;
    e = r*cos(mtheta);
    n = r*sin(mtheta);

    r0index = find(r>=rbox,1,'first');
    xbox = rres*[-xlen/2 -xlen/2 xlen/2-1 xlen/2-1 -xlen/2];                                               
    ybox = [r(r0index) r(r0index+ylen-1) r(r0index+ylen-1) r(r0index) r(r0index)];
    xvec = linspace(min(xbox),max(xbox),xlen);                                                              
    yvec = linspace(min(ybox),max(ybox),ylen);
    [xgrid, ygrid] = meshgrid(xvec,yvec);

    [ebox, nbox] = coordinate_transformation(xbox,ybox,0,0,thetabox.*pi/180);
    [egrid, ngrid] = coordinate_transformation(reshape(xgrid,[],1)',reshape(ygrid,[],1)',0,0,thetabox*pi/180);
    egrid = double(reshape(egrid,size(xgrid)));
    ngrid = double(reshape(ngrid,size(ygrid)));
    [lonbox, latbox] = km2lonlat(lonp(1,1), latp(1,1), ebox./1000, nbox./1000);
    [longrid, latgrid] = km2lonlat(lonp(1,1), latp(1,1), egrid./1000, ngrid./1000);

    if ns==1
        save([framepath 'lonlatgrid.mat'],'longrid','latgrid');
    end
    
    % GPU array setup
    % input grid
    gr = gpuArray(double(r(r1:r2)));
    gtheta = gpuArray(double(unwrap(theta(az1-minf1turning:az2-minf1turning))));
    % output grid
    grgridout = gpuArray(sqrt(egrid.^2+ngrid.^2));
    gthetagridout = gpuArray(pi/2-atan2(ngrid,egrid));
    gthetagridout(gthetagridout<pi/2) = gthetagridout(gthetagridout<pi/2) + 2*pi;
    if ploton
        figure(1)
%         subplot(1,2,1)
        pcolor(e,n,(nanmean(double(f1mag),3)));
        shading flat
        axis square
        colorbar
        caxis([0 10])
        axis([-5000 5000 -5000 5000]);
        title('Frequency 1 Magnitude');
        line([r(r1).*cos(pi/2 - theta(az1-minf1turning)),r(r2).*cos(pi/2 - theta(az1-minf1turning))],[r(r1).*sin(pi/2 - theta(az1-minf1turning)),r(r2).*sin(pi/2 - theta(az1-minf1turning))],'color','r');
        line([r(r1).*cos(pi/2 - theta(az2-minf1turning)),r(r2).*cos(pi/2 - theta(az2-minf1turning))],[r(r1).*sin(pi/2 - theta(az2-minf1turning)),r(r2).*sin(pi/2 - theta(az2-minf1turning))],'color','r');
        line(ebox,nbox,'color','w');
        drawnow;

        
%         subplot(1,2,2)
%         pcolor(x(r1:r2,az1:az2),y(r1:r2,az1:az2),nanmean(double(f1phi(r1:r2,az1:az2,2:2:end))/4096,3));
% %         pcolor(x(:,az1:az2),y(:,az1:az2),nanmean(double(f1time(1:r2,az1:az2,:)),3));
%         shading flat
%         axis square
%         colorbar
%         caxis([-pi pi]);
%         axis([-14000 0 -7000 7000]);
%         title('Frequency 1 Phase');
%         drawnow;
        
        figure(11)
        subplot(1,2,1)
        hist(reshape(f1mag,[],1),1000);
        xlim([0 50]);
        title(['F1 Magnitude - mean = ' num2str(nanmean(reshape(f1mag,[],1)),'%3.2f') ' sigma = ' num2str(nanstd(reshape(f1mag,[],1)),'%3.2f')]);
        subplot(1,2,2)
        hist(reshape(f1phi,[],1)./4095,[-pi:pi/100:pi]);
        xlim([-pi pi]);
        title('F1 Phase');
        drawnow;
        
        figure(2)
%         subplot(1,2,1)
        pcolor(e,n,nanmean(double(f2mag),3));
        shading flat
        axis square
        colorbar
        caxis([0 10])
        axis([-5000 5000 -5000 5000]);
        title('Frequency 2 Magnitude');
        line([r(r1).*cos(pi/2 - theta(az1-minf1turning)),r(r2).*cos(pi/2 - theta(az1-minf1turning))],[r(r1).*sin(pi/2 - theta(az1-minf1turning)),r(r2).*sin(pi/2 - theta(az1-minf1turning))],'color','r');
        line([r(r1).*cos(pi/2 - theta(az2-minf1turning)),r(r2).*cos(pi/2 - theta(az2-minf1turning))],[r(r1).*sin(pi/2 - theta(az2-minf1turning)),r(r2).*sin(pi/2 - theta(az2-minf1turning))],'color','r');
        line(ebox,nbox,'color','w');
        drawnow;
 
%         subplot(1,2,2)
%         pcolor(x(r1:r2,az1:az2),y(r1:r2,az1:az2),nanmean(double(f2phi(r1:r2,az1:az2,2:2:end))/4096,3));
%         shading flat
%         axis square
%         colorbar
%         caxis([-pi pi]);
%         axis([-14000 0 -7000 7000]);
%         title('Frequency 2 Phase');
%         drawnow;
        
        figure(12)
        subplot(1,2,1)
        hist(reshape(f2mag,[],1),4000);
        xlim([0 50]);
        title(['F2 Magnitude - mean = ' num2str(nanmean(reshape(f2mag,[],1)),'%3.2f') ' sigma = ' num2str(nanstd(reshape(f2mag,[],1)),'%3.2f')]);
        subplot(1,2,2)
        hist(reshape(f2phi,[],1)./4095,[-pi:pi/100:pi]);
        xlim([-pi pi]);
        title('F2 Phase');
        drawnow;
%     %     saveas(gcf,[datapath 'scan' d1(ns).name(4:end-4) '.jpg'],'jpg');
    end
    
    % GPU interp2 F1
    gmag = gpuArray(nansum(double(f1mag(r1:r2,az1-minf1turning:az2-minf1turning,:)),3));
    gtime = gpuArray(nanmean(f1time(r1:r2,az1-minf1turning:az2-minf1turning,:),3));

    tic
    gframe = interp2(gtheta,gr,gmag,gthetagridout,grgridout);
    gframetime = interp2(gtheta,gr,gtime,gthetagridout,grgridout);
    f1frame = gather(gframe);
    f1frametime = gather(gframetime);
    toc
    
    if sum(sum(isnan(f1frametime)>0))
        disp('F1 NaN')
        f1frame(isnan(f1frame)) = 0;
    end

    figure(3)
    pcolor(longrid, latgrid, f1frame);
    shading flat
    axis square
    caxis([0 100])
    title('F1')
    drawnow

    disp('saving F1 file');
    save([framepath 'F1_' datestr(scandnum,'yyyymmddHHMMSS') '.mat'],'f1frame','f1frametime','scandnum');
       
    % GPU interp2 F2
    gmag = gpuArray(nansum(double(f2mag(r1:r2,az1-minf1turning:az2-minf1turning,:)),3));
    gtime = gpuArray(nanmean(f2time(r1:r2,az1-minf1turning:az2-minf1turning,:),3));

    tic
    gframe = interp2(gtheta,gr,gmag,gthetagridout,grgridout);
    gframetime = interp2(gtheta,gr,gtime,gthetagridout,grgridout);
    f2frame = gather(gframe);
    f2frametime = gather(gframetime);
    toc
    
    if sum(sum(isnan(f2frametime)>0))
        disp('F2 NaN')
        f2frame(isnan(f2frame)) = 0;
    end

    figure(4)
    pcolor(longrid, latgrid, f2frame);
    shading flat
    axis square
    caxis([0 100])
    title('F2')
    drawnow

    disp('saving F2 file');
    save([framepath 'F2_' datestr(scandnum,'yyyymmddHHMMSS') '.mat'],'f2frame','f2frametime','scandnum');
            
%     R = xcorr2(f1frame,f2frame)./sqrt(dot(reshape(f1frame,[],1),reshape(f1frame,[],1)).*dot(reshape(f2frame,[],1),reshape(f2frame,[],1)));
%     figure
%     pcolor(R);
%     shading flat
%     axis square
%     [v,r,c] = max2d(R);
%     title(['Max = ' num2str(v,'%4.2f') ' at ' num2str(r,'%d') ',' num2str(c,'%d')]); 
end
      