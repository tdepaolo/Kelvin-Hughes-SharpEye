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

for ns=1:length(index)-1
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
        hold off
%         line(ebox,nbox,'color','w');
        
        figure(2)
        pcolor(ep,np,nansum(double(f2mag),3));
%         pcolor(ep,np,double(f2mag(:,:,1)+f2mag(:,:,2)));
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
    % sampling rate (doing frequencies separately)
    fs = 1/(double(fileheader.fsamples)*double(fileheader.sampletime)*1e-9);
    % sample time
    tsamp = (1/fs)*1e6;   % put it in usec

    naz = floor((az2-az1)/squint);
    Dfft1 = NaN.*ones(r2,naz);
    Dfft2 = NaN.*ones(r2,naz);
%     SNR1 = NaN.*ones(r2,naz);
%     SNR2 = NaN.*ones(r2,naz);
    CONF1 = NaN.*ones(r2,naz);
    CONF2 = NaN.*ones(r2,naz);
    tdop1 = NaN.*ones(r2,naz);
    tdop2 = NaN.*ones(r2,naz);
    f1mmag = NaN.*ones(r2,naz);
    f1smag = NaN.*ones(r2,naz);
    f2mmag = NaN.*ones(r2,naz);
    f2smag = NaN.*ones(r2,naz);
    
    % output theta vector
    thetadop = NaN.*ones(1,naz);
    
    azcount = 1;
    
    for nth=az1:squint:az2-squint
        tic
        for nr=r1:r2
            if npulses==1
                % magnitude
                f12m = cat(1,reshape(squeeze(f1mag(nr,nth:nth+squint-1,2:2:end))',[],1),...
                             reshape(squeeze(f2mag(nr,nth:nth+squint-1,2:2:end))',[],1));
                % phase
                f12p = cat(1,reshape(squeeze(f1phi(nr,nth:nth+squint-1,2:2:end))'/4095,[],1),...
                             reshape(squeeze(f2phi(nr,nth:nth+squint-1,2:2:end))'/4095,[],1));
                % time
                f12t = cat(1,reshape(squeeze(f1time(nr,nth:nth+squint-1,2:2:end))',[],1),...
                             reshape(squeeze(f2time(nr,nth:nth+squint-1,2:2:end))',[],1));
            else
                % magnitude
                f1m = reshape(squeeze(f1mag(nr,nth:nth+squint-1,:))',[],1);
                f2m = reshape(squeeze(f2mag(nr,nth:nth+squint-1,:))',[],1);
                % phase
                f1p = reshape(squeeze(f1phi(nr,nth:nth+squint-1,:))'/4095,[],1);
                f2p = reshape(squeeze(f2phi(nr,nth:nth+squint-1,:))'/4095,[],1);
                % time
                f1t = reshape(squeeze(f1time(nr,nth:nth+squint-1,:))',[],1);
                f2t = reshape(squeeze(f2time(nr,nth:nth+squint-1,:))',[],1);
            end

            % remove NaNs 
            idx = ~isnan(f1m);
            f1m = f1m(idx);
            f1p = f1p(idx);
            f1t = f1t(idx);
            idx = ~isnan(f2m);
            f2m = f2m(idx);
            f2p = f2p(idx);
            f2t = f2t(idx);
            
            f1mmag(nr,azcount) = mean(f1m);
            f2mmag(nr,azcount) = mean(f2m);
            f1smag(nr,azcount) = std(f1m);
            f2smag(nr,azcount) = std(f2m);
            
            % Doppler frequency
            % inside of range cell 205, there are no samples for p2 f1, and the f1 arrays may be empty
            if ~isempty(f1m)
                [f1Dfft,f1Dpp,conf1,f2Dfft,f2Dpp,conf2] = compute_Doppler_frequency_separate(f1m,f1p,f1t,f2m,f2p,f2t,tsamp);
                Dfft1(nr,azcount) = f1Dfft;
                Dfft2(nr,azcount) = f2Dfft;
%                 SNR1(nr,azcount) = snr1;
%                 SNR2(nr,azcount) = snr2;
                CONF1(nr,azcount) = conf1;
                CONF2(nr,azcount) = conf2;
                tdop1(nr,azcount) = mean(f1t);
                tdop2(nr,azcount) = mean(f2t);
            end
        end
        toc
        st = sprintf('SharpEye_Doppler_separate: azimuth bin %d', azcount);
        disp(st);
        thetadop(azcount) = mean(theta(nth:nth+squint));
        azcount = azcount + 1;  
    end
    
    % reduce to desired ranges
    Dfft1 = Dfft1(r1:r2,:);
    Dfft2 = Dfft2(r1:r2,:);
%     SNR1 = SNR1(r1:r2,:);
%     SNR2 = SNR2(r1:r2,:);
    CONF1 = CONF1(r1:r2,:);
    CONF2 = CONF2(r1:r2,:);
    tdop1 = tdop1(r1:r2,:);
    tdop2 = tdop2(r1:r2,:);
%     f12tdop = tdop(:,:);
    f1mmag = f1mmag(r1:r2,:);
    f1smag = f1smag(r1:r2,:);
    f2mmag = f2mmag(r1:r2,:);
    f2smag = f2smag(r1:r2,:);
    
%     esub = e(r1:r2,az1:squint:az2);
%     nsub = n(r1:r2,az1:squint:az2);
    
    f1vfft = sol.*Dfft1(:,:)./(2*txf1);
    f2vfft = sol.*Dfft2(:,:)./(2*txf2);
%     f12vpp = sol.*Dpp(:,:)./(2*((txf1+txf2)/2));

    lonsub = lonp(r1:r2,az1:squint:az2-squint);
    latsub = latp(r1:r2,az1:squint:az2-squint);
    
    figure(3)
    subplot(1,2,1)
    pcolor(lonsub,latsub,f1vfft)
    axis([min2d(lonsub) max2d(lonsub) min2d(latsub) max2d(latsub)]);
    caxis([-3 3]);
%     colorbar;
    axis square equal
    shading interp
    title('Doppler Velocity');
    subplot(1,2,2)
    hist(reshape(f1vfft,[],1),(-5:0.1:5));
    xlim([-5 5]);
    title(['F1 FFT Doppler Velocity Distribution - mean = ' num2str(nanmean(reshape(f1vfft,[],1)),'%3.1f') ' sigma = ' num2str(nanstd(reshape(f1vfft,[],1)),'%3.1f')]);
    drawnow
    
    figure(4)
    subplot(1,2,1)
    pcolor(lonsub,latsub,f2vfft)
    axis([min2d(lonsub) max2d(lonsub) min2d(latsub) max2d(latsub)]);
    caxis([-3 3]);
%     colorbar;
    axis square equal
    shading flat
    title('Doppler Velocity');
    subplot(1,2,2)
    hist(reshape(f2vfft,[],1),(-5:0.1:5));
    xlim([-5 5]);
    title(['F2 FFT Doppler Velocity Distribution - mean = ' num2str(nanmean(reshape(f2vfft,[],1)),'%3.1f') ' sigma = ' num2str(nanstd(reshape(f2vfft,[],1)),'%3.1f')]);
    drawnow
    
    f12vfft = (f1vfft+f2vfft)/2;
    f12tdop = (tdop1+tdop2)/2;
%     f12snr = (SNR1+SNR2)/2;
    f12conf = (CONF1+CONF2)/2;
    
    figure(5)
    subplot(1,2,1)
    pcolor(lonsub,latsub,f12vfft)
%     axis([-7500 0 -1500 6000]);
    caxis([-2 2]);
    axis square equal
    shading flat
    title('Mean F1/F2 Doppler Velocity');
    subplot(1,2,2)
    hist(reshape(f12vfft,[],1),(-5:0.1:5));
    xlim([-5 5]);
    title(['Mean F1/F2 FFT Doppler Velocity Distribution - mean = ' num2str(nanmean(reshape(f12vfft,[],1)),'%3.1f') ' sigma = ' num2str(nanstd(reshape(f12vfft,[],1)),'%3.1f')]);
    drawnow
    
    figure(6)
    subplot(1,2,1)
    pcolor(lonsub,latsub,CONF1);
    axis square equal
    shading flat
    caxis([0 1])
    title('CONF1');
    subplot(1,2,2)
    pcolor(lonsub,latsub,CONF2);
    axis square equal
    shading flat
    caxis([0 1])
    title('CONF2');
    
%     figure(7)
%     subplot(1,2,1)
%     pcolor(lonsub,latsub,10*log10(SNR1));
%     axis square equal
%     shading flat
% %     caxis([0 1])
%     title('SNR1');
%     subplot(1,2,2)
%     pcolor(lonsub,latsub,10*log10(SNR2));
%     axis square equal
%     shading flat
% %     caxis([0 1])
%     title('SNR2');
    
    Doppler_stats = [Doppler_stats; scandnum nanmean(reshape(f1vfft,[],1)) nanstd(reshape(f1vfft,[],1)) nanmean(reshape(CONF1,[],1)) nanstd(reshape(CONF1,[],1)) ...
                                             nanmean(reshape(f2vfft,[],1)) nanstd(reshape(f2vfft,[],1)) nanmean(reshape(CONF2,[],1)) nanstd(reshape(CONF2,[],1))];
%     mean_SNR1 = 10*log10(nanmean(reshape(SNR1,[],1)))
%     mean_SNR2 = 10*log10(nanmean(reshape(SNR2,[],1)))
    
    save([scanpath 'f12vfft' datestr(scandnum,'yyyymmddHHMMSS') '.mat'],'f12vfft','lonsub','latsub','f12tdop','thetadop','f12conf');
%     save([scanpath 'f12vpp' datestr(scandnum,'yyyymmddHHMMSS') '.mat'],'f12vpp','lonsub','latsub','f12tdop','thetadop');
    
end
save([scanpath 'Doppler_stats.mat'],'Doppler_stats');

