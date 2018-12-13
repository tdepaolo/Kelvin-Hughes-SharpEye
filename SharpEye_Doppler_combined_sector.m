% Comput SharpEye Doppler velocity from both frequencies F1 and F2
clear all
close all

load('parameters.mat');
% data file
filename = uigetfile([datapath '*.dat']);
% dd = dir([datapath '*.dat']);
% load the file index
load([datapath filename(1:6) '_index.mat']);
% di = dir([datapath '*.mat']);
% load([datapath di(3).name]);
% startscan = find(index(:,1)>=datenum([2017 09 13 23 35 0]));
startscan = 1;
Doppler_stats = [];

for ns=92:length(index)-1
    tic
    [scandnum,...
     f1time,f1mag,f1phi,...
     f2time,f2mag,f2phi,...
     fileheader,...
     minf1turning,maxf1turning,...
     minf2turning,maxf2turning,...
     squint,...
     r,theta,...
     lonp,latp] = read_SharpEye_sector_stitch([datapath filename], index(ns,2));
    toc
    
    % find range and azimuth indexes to Doppler process, defines a wedge, from parameters
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
    
    if ploton

        [ep,np]=lonlat2km(khlongitude,khlatitude,lonp,latp);
        ep = ep*1000;
        np = np*1000;
        
        figure(1)
        pcolor(ep,np,(nansum(double(f1mag),3)));
        shading flat
        axis equal square
        colorbar
        title(['Frequency 1 Magnitude  ' datestr(scandnum)]);
        caxis([0 50])
        drawnow;
        hold on
        line([r(r1).*cos(pi/2 - theta(az1)),r(r2).*cos(pi/2 - theta(az1))],[r(r1).*sin(pi/2 - theta(az1)),r(r2).*sin(pi/2 - theta(az1))],'color','r');
        line([r(r1).*cos(pi/2 - theta(az2)),r(r2).*cos(pi/2 - theta(az2))],[r(r1).*sin(pi/2 - theta(az2)),r(r2).*sin(pi/2 - theta(az2))],'color','r');
%         line(ebox,nbox,'color','w');
        
        figure(2)
        pcolor(ep,np,nansum(double(f2mag),3));
        shading flat
        axis equal square
        colorbar
        caxis([0 50])
        title(['Frequency 2 Magnitude  ' datestr(scandnum)]);
        drawnow;
 
        
%         figure(12)
%         subplot(1,2,1)
%         hist([reshape(f1mag,[],1) reshape(f2mag,[],1)],4000);
%         xlim([0 50]);
%         title(['F2 Magnitude - mean = ' num2str(nanmean(reshape(f2mag,[],1)),'%3.2f') ' sigma = ' num2str(nanstd(reshape(f2mag,[],1)),'%3.2f')]);
%         subplot(1,2,2)
%         hist(reshape(f2phi,[],1)./4095,[-pi:pi/100:pi]);
%         xlim([-pi pi]);
%         title('F2 Phase');
%         drawnow;
%     %     saveas(gcf,[datapath 'scan' d1(ns).name(4:end-4) '.jpg'],'jpg');
    end
    
    
    % Doppler frequency - collect all samples within "squint" turnings,
    % from each frequency.  there are approximately 11 turnings per degree, but it's simpler if 
    % you use the squint offset
    
    % set the number of pulses, P1 and P2 = 2, P2 = 1 (don't use P1)
    npulses = 2;
    % sampling rate
    fs = npulses/(double(fileheader.fsamples)*double(fileheader.sampletime)*1e-9);
    % sample time
    tsamp = (1/fs)*1e6;   % put it in usec

    naz = floor((az2-az1)/squint);
    Dfft = NaN.*ones(r2,naz);
    SNR = NaN.*ones(r2,naz);
    tdop = NaN.*ones(r2,naz);
    mmag = NaN.*ones(r2,naz);
    smag = NaN.*ones(r2,naz);
    % output theta vector
    thetadop = NaN.*ones(1,naz);
    
    azcount = 1;
    
    for nth=az1:squint:az2-squint
        tic
        for nr=r1:r2
            if npulses==1
                % magnitude
                f12m = cat(1,reshape(squeeze(f1mag(nr,nth:nth+squint-1,2:2:end)),[],1),...
                             reshape(squeeze(f2mag(nr,nth:nth+squint-1,2:2:end)),[],1));
                % phase
                f12p = cat(1,reshape(squeeze(f1phi(nr,nth:nth+squint-1,2:2:end))/4095,[],1),...
                             reshape(squeeze(f2phi(nr,nth:nth+squint-1,2:2:end))/4095,[],1));
                % time
                f12t = cat(1,reshape(squeeze(f1time(nr,nth:nth+squint-1,2:2:end)),[],1),...
                             reshape(squeeze(f2time(nr,nth:nth+squint-1,2:2:end)),[],1));
            else
                % magnitude
                f12m = cat(1,reshape(squeeze(f1mag(nr,nth:nth+squint-1,:)),[],1),...
                             reshape(squeeze(f2mag(nr,nth:nth+squint-1,:)),[],1));
                % phase
                f12p = cat(1,reshape(squeeze(f1phi(nr,nth:nth+squint-1,:))/4095,[],1),...
                             reshape(squeeze(f2phi(nr,nth:nth+squint-1,:))/4095,[],1));
                % time
                f12t = cat(1,reshape(squeeze(f1time(nr,nth:nth+squint-1,:)),[],1),...
                             reshape(squeeze(f2time(nr,nth:nth+squint-1,:)),[],1));
            end

            % remove NaNs 
            idx = ~isnan(f12m);
            f12m = f12m(idx);
            f12p = f12p(idx);
            f12t = f12t(idx);
            
            mmag(nr,azcount) = mean(f12m);
            smag(nr,azcount) = std(f12m);
            
            % Doppler frequency
            % inside of range cell 205, there are no samples for p2 f1, and the f1 arrays may be empty
            if ~isempty(f12m)
                [f12Dfft, snr] = compute_Doppler_frequency(f12m,f12p,f12t,tsamp);
                Dfft(nr,azcount) = f12Dfft;
                SNR(nr,azcount) = snr;
                tdop(nr,azcount) = mean(f12t);
            end
        end
        toc
        st = sprintf('SharpEye_Doppler_combined: azimuth bin %d', azcount);
        disp(st);
        thetadop(azcount) = mean(theta(nth:nth+squint));
        azcount = azcount + 1;  
    end
    
    % reduce to desired ranges
    Dfft = Dfft(r1:r2,:);
    SNR = SNR(r1:r2,:);
    tdop = tdop(r1:r2,:);
    f12tdop = tdop(:,:);
    mmag = mmag(r1:r2,:);
    smag = smag(r1:r2,:);
    
%     esub = e(r1:r2,az1:squint:az2);
%     nsub = n(r1:r2,az1:squint:az2);
    
    f12vfft = sol.*Dfft(:,:)./(2*((txf1+txf2)/2));
%     f12vpp = sol.*Dpp(:,:)./(2*((txf1+txf2)/2));

    lonsub = lonp(r1:r2,az1:squint:az2-squint);
    latsub = latp(r1:r2,az1:squint:az2-squint);
    
    figure(3)
    subplot(1,2,1)
    pcolor(lonsub,latsub,f12vfft)
    axis([min2d(lonsub) max2d(lonsub) min2d(latsub) max2d(latsub)]);
    caxis([-2 2]);
    colorbar;
    axis square equal
    shading flat
    title('Doppler Velocity');
    subplot(1,2,2)
    hist(reshape(f12vfft,[],1),(-5:0.1:5));
    xlim([-5 5]);
    title(['F1/F2 FFT Doppler Velocity Distribution - mean = ' num2str(nanmean(reshape(f12vfft,[],1)),'%3.1f') ' sigma = ' num2str(nanstd(reshape(f12vfft,[],1)),'%3.1f')]);
    drawnow
    
%     figure(5)
%     subplot(1,2,1)
%     pcolor(lonsub,latsub,f12vpp)
% %     axis([-7500 0 -1500 6000]);
% %     caxis([-1 1]);
%     axis square equal
%     shading flat
%     subplot(1,2,2)
%     hist(reshape(f12vpp,[],1),(-5:0.1:5));
%     xlim([-5 5]);
%     title(['F1/F2 PP Doppler Velocity Distribution - mean = ' num2str(nanmean(reshape(f12vpp,[],1)),'%3.1f') ' sigma = ' num2str(nanstd(reshape(f12vpp,[],1)),'%3.1f')]);
%     drawnow
    
    Doppler_stats = [Doppler_stats; scandnum nanmean(reshape(f12vfft,[],1)) nanstd(reshape(f12vfft,[],1)) 10*log10(nanmean(reshape(SNR,[],1))) 10*log10(nanstd(reshape(SNR,[],1)))];
    mean_SNR = 10*log10(nanmean(reshape(SNR,[],1)))
    
    save([scanpath 'f12vfft' datestr(scandnum,'yyyymmddHHMMSS') '.mat'],'f12vfft','lonsub','latsub','f12tdop','thetadop');
%     save([scanpath 'f12vpp' datestr(scandnum,'yyyymmddHHMMSS') '.mat'],'f12vpp','lonsub','latsub','f12tdop','thetadop');
    
end
