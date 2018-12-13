% generate ocean surface map from Kelvin Hughes radar scans
% the hard part is unfolding the aliased spectrum to set up the dispersion
% filter
% none of this uses the GPU yet...need to work on that....

clear
close all

% get processing parameters
load('parameters.mat');
% df1 = dir([framepath 'F1*.mat']);
% df2 = dir([framepath 'F2*.mat']);
dv = dir([datapath 'f12vfft*.mat']);
radar_time = zeros(length(dv),1);
for nf=1:length(dv)
    radar_time(nf) = datenum(dv(nf).name(8:21),'yyyymmddHHMMSS');
end
% lidar swath times
t1 = datenum([2017 9 20 0 21 04]);
t2 = datenum([2017 9 20 0 23 10]);
tidx1 = find(radar_time >= t1,1,'first');
tidx2 = find(radar_time >= t2,1,'first');

% load([framepath 'lonlatgrid.mat']);

% set up inverted spatial axes, these are wavenumbers (rad/m)
ks = 2*pi/rres;
kxres = ks/(xlen*2);
kyres = ks/(ylen*2);
kx = (kxres*[-(xlen*2)/2:((xlen*2)/2)-1]);
ky = (kyres*[-(ylen*2)/2:((ylen*2)/2)-1]);
[kxout, kyout] = meshgrid(kx,ky);
kmag1 = kx((xlen*2)/2+1:end);

% These must be edited depending on how the spectrum looks
% kcutp1 is the lower cutoff, and kcutp2 is the upper cutoff that define
% the wavenumber passband for the shifted dispersion filter
kcutp1 = find(kmag1 >= kcutoff_low, 1 );
kcutp2 = find(kmag1 >= kcutoff_hi, 1);

% angle vector
tres = 2*pi/(xlen*2);
theta = (tres*[0:(xlen*2)-1]);
% polar coordinates
xp = kmag1'*cos(theta);
yp = kmag1'*sin(theta);
% dispersion curve as a funtion of w
dispersionw = sqrt(g.*kmag1.*tanh(kmag1*h));
% dispersionw2 = 2*sqrt(g./2*kmag1.*tanh(kmag1*h/2));

% setup arrays for shifted dispersion curve fit
A=[dispersionw(kcutp1:kcutp2).' kmag1(kcutp1:kcutp2).'];
kindex = kcutp1:kcutp2;
inkmag = kmag1(kindex);
b = zeros(length(kindex),1);
smax = zeros(length(kindex),1);

% 2D wavenumber array and dispersion curve
kmag2 = sqrt(kxout.^2 + kyout.^2);
dispersion2w = sqrt(g.*kmag2.*tanh(kmag2*h));
% dispersion2w2 = 2*sqrt(g./2*kmag2.*tanh(kmag2*h/2));
   
mask = zeros((xlen*2),(ylen*2),(tlen*3));
c = zeros(xlen*2,ylen*2,tlen*2);
t = zeros(xlen,ylen,tlen);
cdnum = zeros(tlen,1);

startframe = 1;
endframe = length(dv);

for nframe=startframe:endframe
    tic
    
    load([datapath dv(nframe).name]);
    [e,n] = lonlat2km(khlongitude,khlatitude,lonsub,latsub);
    e = e.*1000;
    n = n.*1000;
%     
%     figure(1)
%     pcolor(lonsub,latsub,f12vfft)
%     shading interp
%     axis equal
%     caxis([-2 2]);
%     hold on
%     load('SPOT-0026_gps.mat')
%     plot(lon,lat,'.r');
%     figure(2)
%     pcolor(e,n,f12vfft)
%     shading interp
%     axis equal
%     caxis([-2 2]);
    
    if nframe==startframe
        
        r = rres.*(0:9999)';
        r1 = find(r>=rmin,1,'first');
        r2 = find(r>=rmax,1,'first');
        r = r(r1:r2);
%         e = r*cos(pi/2-thetadop);
%         n = r*sin(pi/2-thetadop);
        
        r0index = find(r>=rbox,1,'first');
        xbox = rres*[-xlen/2 -xlen/2 xlen/2-1 xlen/2-1 -xlen/2];                                               
        ybox = [r(r0index) r(r0index+ylen-1) r(r0index+ylen-1) r(r0index) r(r0index)];
        xvec = linspace(min(xbox),max(xbox),xlen);                                                              
        yvec = linspace(min(ybox),max(ybox),ylen);
        [xgrid, ygrid] = meshgrid(xvec,yvec);

        [ebox, nbox] = coordinate_transformation(xbox,ybox,0,0,thetabox.*pi/180);
        [egrid, ngrid] = coordinate_transformation(reshape(xgrid,[],1)',reshape(ygrid,[],1)',0,0,thetabox*pi/180);

        % rotate again to desired orientation
        [ebox, nbox] = coordinate_rotation(ebox,nbox,r(r0index)*cos(pi/2 - thetabox*pi/180),r(r0index)*sin(pi/2 - thetabox*pi/180),dthetabox*pi/180);
        [rxgrid, rygrid] = coordinate_rotation(egrid,ngrid,r(r0index)*cos(pi/2 - thetabox*pi/180),r(r0index)*sin(pi/2 - thetabox*pi/180),dthetabox*pi/180);
        egrid = reshape(rxgrid,size(xgrid));
        ngrid = reshape(rygrid,size(ygrid));
        
        % compute angles to the grid points (Cartesian, not compass angles)
%         grangles = atan2(ngrid,egrid);
        grangles = mod(atan2(ngrid,egrid),2*pi);

        [lonbox, latbox] = km2lonlat(khlongitude, khlatitude, ebox./1000, nbox./1000);
        [longrid, latgrid] = km2lonlat(khlongitude, khlatitude, egrid./1000, ngrid./1000);
        
%         hold on
%         plot(ebox,nbox,'w');
%         hold off
%         drawnow;
        save([wavepath 'box.mat'],'lonbox','latbox');
        save([wavepath 'lonlatgrid.mat'],'longrid','latgrid');
        
        % GPU array setup
        % input grid
        gr = gpuArray(r);
        gtheta = gpuArray(thetadop);

        % output grid
        grgridout = gpuArray(sqrt(egrid.^2+ngrid.^2));
        gthetagridout = gpuArray(pi/2-atan2(ngrid,egrid));
        gthetagridout(gthetagridout<pi/2) = gthetagridout(gthetagridout<pi/2) + 2*pi;

    end

    figure(1)
%     subplot(1,2,1)
    pcolor(e,n,f12vfft);
    shading interp
    axis equal
    caxis([-2 2])
    hold on
    plot(ebox,nbox,'k');
    hold off
    colorbar
%     subplot(1,2,2)
%     pcolor(e,n,f12conf);
%     shading interp
%     axis equal
%     caxis([0 1]);
%     colorbar
    drawnow;
    
    % GPU interp2
    gv = gpuArray(f12vfft);
    gtime = gpuArray(f12tdop);

    tic
    gframe = interp2(gtheta,gr,gv,gthetagridout,grgridout);
    gframetime = interp2(gtheta,gr,gtime,gthetagridout,grgridout);
    frame = gather(gframe);
    frametime = gather(gframetime);
    toc
    
%     figure(2)
%     pcolor(egrid,ngrid,frame);
%     shading flat
%     axis square tight
%     caxis([-2 2])
%     drawnow;
    
    % scale the radial velocity according to the relative angle of the
    % dominant wave angle, where the box is pointing
    % convert to Cartesian angle
    wangle = mod(90-(thetabox+dthetabox),360)*pi/180;
%     sframe = frame./cos(wangle-grangles);
    sframe = frame;
    
    figure(2)
    pcolor(egrid,ngrid,sframe);
    shading flat
    axis square tight
    caxis([-2 2])
    drawnow;
    
    if (nframe<=startframe+tlen-1)
        % stack up tlen frames
        c(xlen/2+1:xlen/2+xlen,ylen/2+1:ylen/2+ylen,tlen/2+nframe-startframe+1) = sframe;
        t(:,:,nframe-startframe+1) = frametime;
        cdnum(nframe-startframe+1) = datenum(dv(nframe).name(8:21),'yyyymmddHHMMSS');      
    else
        % put the next frame on the bottom 
        c(xlen/2+1:xlen/2+xlen,ylen/2+1:ylen/2+ylen,tlen/2+1:tlen/2+tlen-1) =....
            c(xlen/2+1:xlen/2+xlen,ylen/2+1:ylen/2+ylen,tlen/2+2:tlen/2+tlen);
        c(xlen/2+1:xlen/2+xlen,ylen/2+1:ylen/2+ylen,tlen/2+tlen) = sframe;
        t(:,:,1:end-1) = t(:,:,2:end);
        t(:,:,end) = frametime;
        cdnum(1:end-1) = cdnum(2:end);
        cdnum(end) = datenum(dv(nframe).name(8:21),'yyyymmddHHMMSS');
    end

    if nframe >= startframe+tlen-1
                
        % 3-D fft
        % determined the average rotation (time sample) rate using one
        % pixel location
        Tsampvec = datevec(mean(diff(cdnum + datenum(squeeze(t(xlen,ylen,:))/(24*60*60*1e6)))));
        Tsamp = Tsampvec(6);
        
        % set up inverted temporal axes, these are in radian frequency (rad/s)
        ws = 2*pi/Tsamp;
        wres = ws/(tlen*2);
        w = (wres*[-(tlen*2)/2:((tlen*2)/2)-1]);
        pw = w((tlen*2)/2+1:end);    
        % make triple frequency spectrum for de-aliasing
        pw2 = cat(2,pw,pw(end)+(wres.*(1:tlen)));
        pw3 = cat(2,pw2,pw2(end)+(wres.*(1:tlen)));
        w3 = (wres*[-(tlen*3):tlen*3-1]);
        % low frequency cutoff index
        wcutp = find(pw3 > wcutoff, 1 );
        
        disp('FFT');
        c(isnan(c)) = 0;
        dspectrum = fftshift(fftn(c));
        
%         figure(19)
%         for nn=1:size(dspectrum,3)
%             pcolor(kxout,kyout,10*log10(abs(dspectrum(:,:,nn))));
%             shading flat
%             axis([-0.2 0.2 -0.2 0.2]);
%             title(['w = ' num2str(w(nn),'%4.3f')]);
%             drawnow
%         end
            
        % make triple frequency spectrum, postitive frequencies only
        % negative spectrum is complex conjugate of positive, and includes
        % and extra frequency bin
        % due to aliasing, positive wavenumber energy shows up at aliased
        % negative frequencies, so to make the spectrum, add the spectral
        % slice from each frequency to itself, flipped up/down and
        % left/right, and average the two
        pdspectrum = zeros(size(dspectrum,1),size(dspectrum,2),3*size(dspectrum,3)./2);
        pdspectrum(:,:,1:tlen) = dspectrum(:,:,tlen+1:tlen*2);
        pdspectrum(:,:,tlen+1:tlen*2) = dspectrum(:,:,1:tlen);
        pdspectrum(:,:,tlen*2+1:tlen*3) = dspectrum(:,:,tlen+1:tlen*2);
        
        
%         figure(19)
%         for nn=1:size(pdspectrum,3)
%             pcolor(kxout,kyout,10*log10(abs(pdspectrum(:,:,nn))));
%             shading flat
%             axis([-0.2 0.2 -0.2 0.2]);
%             title(['w = ' num2str(pw3(nn),'%4.3f')]);
%             caxis([0 50]);
%             colorbar
%             drawnow
%         end
        % compute power spectrum for SNR and detection threshold 
        psd2 = 2*pdspectrum.*conj(pdspectrum);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Determine Doppler shifted dispersion relationship, due to current
        % for each positive frequency, change (kx,ky,w) coordinates to polar (|k|,theta,w)
        % integrate over theta to get a two dimensional power spectral density S(|k|,w)
        disp('Doppler Shift');
   
        ippsd = zeros(length(pw3),length(kmag1));
        
        % converts 3D PSD to polar coordinates, integrates out
        % theta to make a 2D spectrum
        for nn=1:length(pw3)
            ppsd = interp2(kxout,kyout,psd2(:,:,nn),xp,yp,'linear',max(kx));
            ippsd(nn,:) = sum(abs(ppsd),2)/size(ppsd,2)';
        end
        
        % plot k-w spectrum
        if ploton
            figure(3)
            contourf(kmag1,pw3,10*log10(abs(ippsd)),30);
%             surf(kmag1,pw3,10*log10(abs(ippsd3)));
            axis([0 max(kx) 0 max(pw3)]);
            xlabel('|k| (rad/m)');
            ylabel('omega (rad/s)');
            set(gca,'YTick',pw3);
%             set(gca,'XTick',kmag1);
%             set(gca,'YTick',[0 pw(length(pw)/4) pw(length(pw)/2) pw(3*length(pw)/4) pw(end)]);
            colorbar;
            xlim([0 0.2]);
%             caxis([0 70]);
            hold on;
            line([kmag1(1) kmag1(end)], [pw3(tlen) pw3(tlen)], 'Color', 'r')
            line([kmag1(1) kmag1(end)], [pw3(tlen+1) pw3(tlen+1)], 'Color', 'r')
            line([kmag1(1) kmag1(end)], [pw3(tlen+2) pw3(tlen+2)], 'Color', 'r')
            
            plot(kmag1,dispersionw,'-k');

            plot(kmag1,dispersionw -  4*wres,'-k');
            plot(kmag1,dispersionw +  4*wres,'-k');

            line([kmag1(kcutp1) kmag1(kcutp1)], [pw(1) pw(end)], 'Color', 'k');
            line([kmag1(kcutp2) kmag1(kcutp2)], [pw(1) pw(end)], 'Color', 'k');
            line([kmag1(1) kmag1(end)], [pw(wcutp) pw(wcutp)], 'Color', 'k');
        end
        
        % this reduces the search area on the spectrum within +/- the expected
        % maximum current parameter
        dmask = zeros(size(ippsd));
        for iw = wcutp:length(pw3)
%                 dmask(iw,:) = (pw3(iw) >= dispersionw - (kmag1.*cmax) & ...
%                                (pw3(iw) <= dispersionw + (kmag1.*cmax)));
                dmask(iw,:) = (pw3(iw) >= dispersionw - 4*wres & ...
                               (pw3(iw) <= dispersionw + 4*wres));
        end
        dippsd = ippsd.*dmask;
        dippsd(dippsd==0) = NaN;
        
        % erase part of  the aliased spectrum before filtering
        % this part needs to be adjusted depending on how the spectrum
        % looks
        
        % first aliasing point
        ka1 = pw3(tlen)^2./g;
        ka1idx = find(kmag1>=ka1,1);
        % mark the first aliasing point
        if ploton
            hold on
            plot(kmag1(ka1idx),pw3(tlen),'or')
        end
        % second aliasing point
        ka2 = pw3(2*tlen)^2./g;
        ka2idx = find(kmag1>=ka2,1);
        % mark the second aliasing point
        if ploton
            plot(kmag1(ka2idx),pw3(tlen*2),'or')
            hold off
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % erase unusable spectrum
        % low frequency, high wavenumber, ENABLE/ADJUST THIS
%         dippsd(1:tlen,ka1idx+1:end)=NaN;

        % mid-frequency band, all wavenumbers ENABLE/ADJUST THIS
        dippsd(2*tlen+1:end,1:ka2idx-1) = NaN;
        dippsd(1:2*tlen,ka2idx+1:end) = NaN;
        
        % low-mid frequency, high wavenumber, ENABLE/ADJUST THIS
%         dippsd(1:tlen*2,ka2idx+2:end) = NaN;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
        % plot unfolded spectrum
        if ploton
            figure(4)
            contourf(kmag1,pw3,10*log10(abs(dippsd)),60);
            axis([0 max(kx) 0 max(pw3)]);
            xlabel('|k| (rad/m)');
            ylabel('omega (rad/s)');
    %             set(gca,'YTick',pw);
    %             set(gca,'XTick',kmag1);
%             set(gca,'YTick',[0 pw(length(pw)/4) pw(length(pw)/2) pw(3*length(pw)/4) pw(end)]);
            colorbar;
            xlim([0 0.2]);
%             caxis([0 20]);
            hold on;
            plot(kmag1,dispersionw,'-k');
            plot(kmag1,dispersionw - 4*wres,'-k');
            plot(kmag1,dispersionw + 4*wres,'-k');
        end
       
        % plot red asterisks that will be used for the curve fit, make sure
        % they are on/near the shifted dispersion energy and not in the
        % noise or aliased spectrum
        
        [wmax,iwmax] = max(dippsd(wcutp:end,kindex));
        iwmax = iwmax+wcutp-1;
        if ploton
            hold on
            plot(kmag1(kindex),pw3(iwmax),'r*');
            hold off
        end 
        b = pw3(iwmax).';
        smax = wmax;
        
        inkmag = kmag1(kindex);
        % least squares curve fit for shifted dispersion curve
        x = lscov(A,b);
        % the red line should go through the red asterisks
        if ploton
            hold on
            plot(inkmag, A*x,'-r');
            hold off
        end
      
        % average signal region of S(|k|,w) to find average signal+noise power
        % spectral density
        Ps = mean(smax);
        [Psmax,sidx] = max(smax);
        % average noise region of S(|k|,w) to find average noise power spectral
        % density
        Pn = mean(mean(ippsd(wcutp:wcutp+(tlen*2)/4,length(kmag1)/2:end)));
        if ploton
            figure(3)
            hold on
            line([kmag1(length(kmag1)/2), kmag1(length(kmag1)/2), kmag1(end), kmag1(end), kmag1(length(kmag1)/2)],...
                [pw3(wcutp), pw3(wcutp+(tlen*2)/4), pw3(wcutp+(tlen*2)/4), pw3(wcutp), pw3(wcutp)],...
                'Color','k');
            hold off
        end

        % average SNR
        SNR = ((Ps - Pn)/Pn);
        % maximum SNR
        SNRmax = ((Psmax - Pn)/Pn)
        if ploton
            title(['Maximum SNR ' num2str(10*log10(SNRmax)) ' dB']); 
            drawnow;
        end
        
        % SET A BREAKPOINT HERE TO EXAMINE THE SPECTRUM AND ADJUST
        % PARAMETERS
    
        disp('Filter');
        % 3-D brick wall band pass filter
        % around the dispersion relationship, shifted for current,
        % set all the rest = 0
        snr = (psd2 - Pn)./Pn;
        threshold = SNRmax*0.1; 
                
        Ax = reshape([reshape(dispersion2w,[],1) reshape(kmag2,[],1)]*x, (xlen*2), (ylen*2));
        % upper positive w limit, dispersion+current
        upwdc = Ax + 1*wres;
        % lower positive w limit, dispersion-current
        lpwdc = Ax - 1*wres;
        
        kxin = [];
        kyin = [];
        win = [];
        
        % PARFOR, this does the current estimate, which uses a very narrow
        % band pass filter (mask)
        %                           kxout>=-0.08 & kxout<=0.08 &...
        %                           kyout>0 &...
        for iw=wcutp:length(pw3)
            tmpmask= (lpwdc <= pw3(iw) &...
                          pw3(iw) <= upwdc &...
                          sqrt(kxout.^2+kyout.^2) >= min(inkmag) &...
                          sqrt(kxout.^2+kyout.^2) <= max(inkmag) &...
                          kxout>=-0.08 & kxout<=0.08 &...
                          kyout>0 &...
                          snr(:,:,iw) >= threshold);
                      
            mask(:,:,iw) = gather(tmpmask);

            % use these points for current regression
            if sum(sum(mask(:,:,iw))) ~= 0
                kxin = [kxin; kxout(logical(mask(:,:,iw)))];
                kyin = [kyin; kyout(logical(mask(:,:,iw)))];
                win = [win; pw3(iw).*ones(sum(sum(mask(:,:,iw))),1)];
            end
        end

        % Current vector regression
        % don't want all the wave (kx,ky,w) triplets for current
        % estimation, only the ones on the shifted dispersion line (A*x)
            
        [u,uint,r] = regress(win-sqrt(g*sqrt((kxin.^2+kyin.^2)).*tanh(sqrt(kxin.^2+kyin.^2).*h)),[kxin kyin]);

        % plot the current result
        kxfit = min(kxin):0.001:max(kxin);
        kyfit = min(kyin):0.001:max(kyin);
        [KXFIT, KYFIT] = meshgrid(kxfit,kyfit);
        WFIT = sqrt(g*sqrt(KXFIT.^2+KYFIT.^2).*tanh(sqrt(KXFIT.^2+KYFIT.^2).*h)) + u(1).*KXFIT + u(2).*KYFIT;
        if ploton
            figure(20)
            plot3(kxin,kyin,win,'.');
            hold on
            mesh(KXFIT,KYFIT,WFIT);
            hold off
            title(['Current Estimate: Ux = ' num2str(u(1),'%4.2f') '  Uy = ' num2str(u(2),'%4.2f')]);
            drawnow
        end 
       
        % try to minimize functional
        fun = @(U)sum((win-sqrt(g*sqrt((kxin.^2+kyin.^2)).*tanh(sqrt(kxin.^2+kyin.^2).*h))-(kxin*U(1)+kyin*U(2))).^2);
        U0 = [0,0];
        U = fminsearch(fun,U0)
        
        % now do it again, opening up filter width in frequency for waves
        Ax = reshape([reshape(dispersion2w,[],1) reshape(kmag2,[],1)]*x, (xlen*2), (ylen*2));
        % upper positive w limit, dispersion+current
        upwdc = Ax + 2*wres;
        % lower positive w limit, dispersion-current
        lpwdc = Ax - 2*wres;
            
        kxin = [];
        kyin = [];
        win = [];
        
        % PARFOR, wider bandpass filter (mask)
        %                           kxout>=-0.08 & kxout<=0.08 &...
        %                           kyout>0 &...
        for iw=wcutp:length(pw3)
            tmpmask= (lpwdc <= pw3(iw) &...
                          pw3(iw) <= upwdc &...
                          sqrt(kxout.^2+kyout.^2) >= min(inkmag) &...
                          sqrt(kxout.^2+kyout.^2) <= max(inkmag) &...
                          kxout>=-0.08 & kxout<=0.08 &...
                          kyout>0.001 &...
                          snr(:,:,iw) >= threshold);
                      
            mask(:,:,iw) = gather(tmpmask);

            % use these points for wave inversion
            if sum(sum(mask(:,:,iw))) ~= 0
                kxin = [kxin; kxout(logical(mask(:,:,iw)))];
                kyin = [kyin; kyout(logical(mask(:,:,iw)))];
                win = [win; pw3(iw).*ones(sum(sum(mask(:,:,iw))),1)];
            end
        end
        
        if ploton
            kxfit = min(kxin):0.001:max(kxin);
            kyfit = min(kyin):0.001:max(kyin);
            [KXFIT, KYFIT] = meshgrid(kxfit,kyfit);
            WFIT = sqrt(g*sqrt(KXFIT.^2+KYFIT.^2).*tanh(sqrt(KXFIT.^2+KYFIT.^2).*h)) + u(1).*KXFIT + u(2).*KYFIT; 
            figure(6)
            plot3(kxin,kyin,win,'.');
            hold on
            mesh(KXFIT,KYFIT,WFIT);
            hold off
            axis([min2d(kxout) max2d(kxout) min2d(kyout) max2d(kyout) pw3(1) pw3(end)]);
            xlabel('kx');
            ylabel('ky');
            grid on
            title('Filtered wave components');
            drawnow
        end
        
        pdspectrum = pdspectrum .* mask;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % plot filtered spectrum
        % for each positive frequency, change (kx,ky,w) coordinates to polar (|k|,theta,w)
        % integrate over theta to get a two dimensional power spectral density S(|k|,w)
        psd2 = pdspectrum.*conj(pdspectrum); 
        disp('Filtered Spectrum');
        
        % convert to polar, integrate out theta to make 2D spectrum
        for nn=1:length(pw3)
            ppsd = interp2(kxout,kyout,psd2(:,:,nn),xp,yp,'linear',max(kx));
            ippsd(nn,:) = sum(abs(ppsd),2)/size(ppsd,2)';
        end
        ippsd(ippsd==0)=NaN;

        if ploton
            figure(13)
            contour(kmag1,pw3,10*log10(abs(ippsd)),60);
            axis([0 0.2 0 max(pw3)]);
            xlabel('|k| (rad/m)');
            ylabel('w (rad/s)');
        %     set(gca,'YTick',pw);
            set(gca,'YTick',[0 pw(length(pw)/4) pw(length(pw)/2) pw(3*length(pw)/4) pw(end)]);
            colorbar;

            hold on;
            [row,col]=find(~isnan(ippsd));
            ikmin = min(col);
            ikmax = max(col);
            line([kmag1(ikmin) kmag1(ikmin)], [pw(1) pw(end)]);
            line([kmag1(ikmax) kmag1(ikmax)], [pw(1) pw(end)]);
            hold off;
            drawnow
        end
 
        % take the current out, don't currently use this...just checking to
        % see that (kxin,kyin,wun) is on the dispersion relationship, and
        % it is! 
        wun = win - [kxin kyin]*[u(1); u(2)];

        % apply transfer function from velocity spectrum -> displacment
        % spectrum  (w-dot(k,U))./tanh(|k|*h)
        kangles = mod(90-atan2(kyin,kxin)*180/pi,360);
        tf = (wun./tanh(sqrt(kxin.^2+kyin.^2).*h)).^2;
%         tf = (wun.*cos(kangles*pi/180)).^2;
        pdspectrum(logical(mask)) = pdspectrum(logical(mask))./tf;
        
        disp('IFFT');
        % create conjugate transpose for negative spectrum, concatenate
        % for negative spectrum, eliminate DC (w=0), add a layer of zeros
        % for the most negative frequency
        ndspectrum = conj(flipdim(flipdim(flipdim(pdspectrum(:,:,2:end),1),2),3));
        ndspectrum = cat(3,zeros(xlen*2,ylen*2),ndspectrum);
        % add in the positive spectrum
        dspectrum = cat(3,ndspectrum,pdspectrum(:,:,1:end));
        
%         figure(19)
%         for nn=1:size(dspectrum,3)
%             pcolor(kxout,kyout,10*log10(abs(dspectrum(:,:,nn))));
%             shading flat
%             axis([-0.2 0.2 -0.2 0.2]);
%             title(['w = ' num2str(w3(nn),'%4.3f')]);
%             drawnow
%         end

        zout = ifftn(ifftshift(dspectrum),'symmetric');
        % throw away the zero pad frames
        zout = 2*real(zout(xlen/2+1:xlen/2+xlen,ylen/2+1:ylen/2+ylen,3*tlen/2+1:size(zout,3)-3*tlen/2));
        % the remaining frames have triple the sampling rate
        dnumout = (cdnum(1):datenum((Ts/3)/(24*60*60)):cdnum(end)+datenum((2*Ts/3)/(24*60*60)))';
        tout = zeros(size(zout));
        for nn=1:size(tout,3)
             % fill in the times
            if ~logical(mod(nn-1,3))
                tout(:,:,nn) = t(:,:,(nn-1)/3+1);
                tout(:,:,nn+1) = tout(:,:,nn) + Ts*1e6/3;
                tout(:,:,nn+2) = tout(:,:,nn+1) + Ts*1e6/3;
            end
        end
        
%         figure(6)
%         for nn=1:size(zout,3)
%             if ~logical(mod(nn-1,3))
%                 subplot(1,2,1)
%                 pcolor(egrid,ngrid,c(xlen/2+1:xlen/2+xlen,ylen/2+1:ylen/2+ylen,(tlen/2+1)+(nn-1)/3));
%                 axis equal square;
%                 colormap winter
%                 colorbar;
%                 caxis([-3 3]);
%                 shading interp;
%                 axis tight
%                 title(datestr(cdnum((nn-1)/3+1)));
%             end
%             subplot(1,2,2)
%             pcolor(egrid,ngrid,double(zout(:,:,nn)));
%             axis equal square;
%             colormap winter
%             colorbar;
%             caxis([-Hs Hs]);
%             shading interp;
%             axis tight
%             title(datestr(dnumout(nn)));
%             drawnow
%         end

        % output the central slice
        zout = zout(:,:,3*tlen/2+1);
        tout = tout(:,:,3*tlen/2+1);
        dnumout = dnumout(3*tlen/2+1);

        figure(10)
        pcolor(longrid,latgrid,double(zout));
        axis equal square;
        colormap winter
        colorbar;
        caxis([-Hs Hs]);
        shading interp;
        axis tight
        title(datestr(dnumout));
        drawnow

        fdate = datestr(dnumout,'yyyymmddHHMMSS');

%         save([wavepath 'fs' fdate '.mat'],'kxout','kyout','w','dspectrum','mpassband', 'A', 'b', 'Tsamp' );
        % save the (kx,ky,w) triplets for Doppler velocity inversion
        save([wavepath 'v_wkxky' fdate '.mat'], 'win','kxin','kyin','u');
        % save the waves
        save([wavepath 'v_waves' fdate '.mat'],'zout','tout','dnumout');

    end
end

