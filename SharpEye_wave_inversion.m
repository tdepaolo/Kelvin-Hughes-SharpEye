% generate ocean surface map from Kelvin Hughes radar scans
% the hard part is unfolding the aliased spectrum to set up the dispersion
% filter
% none of this uses the GPU yet...need to work on that....

clear
close all

% get processing parameters
load('parameters.mat');
df1 = dir([framepath 'F1*.mat']);
df2 = dir([framepath 'F2*.mat']);
load([framepath 'lonlatgrid.mat']);

% set up inverted spatial axes, these are wavenumbers (rad/m)
ks = 2*pi/rres;
kxres = ks/(xlen*2);
kyres = ks/(ylen*2);
kx = (kxres*[-(xlen*2)/2:((xlen*2)/2)-1]);
ky = (kyres*[-(ylen*2)/2:((ylen*2)/2)-1]);
[kxout, kyout] = meshgrid(kx,ky);
kmag1 = kx((xlen*2)/2+1:end);

% These must be edited depending on how the spectrum looks
% kcutp1 is the lower cutoff index, and kcutp2 is the upper cutoff index that define
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
b = zeros(length(kindex),1);


% 2D wavenumber array and dispersion curve
kmag2 = sqrt(kxout.^2 + kyout.^2);
dispersion2w = sqrt(g.*kmag2.*tanh(kmag2*h));
% dispersion2w2 = 2*sqrt(g./2*kmag2.*tanh(kmag2*h/2));

mask = zeros((xlen*2),(ylen*2),(tlen*3));
c = zeros(xlen*2,ylen*2,tlen*2);
t = zeros(xlen,ylen,tlen);
cdnum = zeros(tlen,1);
startframe = tlen;
TotFrame = startframe+10;
U = zeros(2*(TotFrame-startframe+1),3);

% just processing one stack to generate (kx,ky,w) triplets
for nframe=1:length(df1)
    tic
    D1 = load([framepath df1(nframe).name]);
    D2 = load([framepath df2(nframe).name]);
    if (nframe<=startframe)
        % it was necessary to output "load" into D if we want to make above
        % for loop a parfoor loop...
        % averaging both frames (from F1 and F2)
        % stack up tlen frames
        c(xlen/2+1:xlen/2+xlen,ylen/2+1:ylen/2+ylen,tlen/2+nframe) = (D2.f2frame + D2.f2frame);
        t(:,:,nframe) = D1.f1frametime;
        cdnum(nframe) = D1.scandnum;      
    else
        % put the next frame on the bottom 
        c(xlen/2+1:xlen/2+xlen,ylen/2+1:ylen/2+ylen,tlen/2+1:tlen/2+tlen-1) =....
            c(xlen/2+1:xlen/2+xlen,ylen/2+1:ylen/2+ylen,tlen/2+2:tlen/2+tlen);
        c(xlen/2+1:xlen/2+xlen,ylen/2+1:ylen/2+ylen,tlen/2+tlen) = (D1.f1frame + D2.f2frame);
        t(:,:,1:end-1) = t(:,:,2:end);
        t(:,:,end) = D1.f1frametime;
        cdnum(1:end-1) = cdnum(2:end);
        cdnum(end) = D2.scandnum;

    end

    if nframe >= startframe
                
        % 3-D fft
%         % zero pad the cube with time slices in the front and back
%         c = cat(3,zeros(xlen,ylen,tlen/2),c);
%         c = cat(3,c,zeros(xlen,ylen,tlen/2));
%         % also zero pad the cube around the spatial area
%         c = cat(1,zeros(xlen/2,ylen,tlen*2),c);
%         c = cat(1,c,zeros(xlen/2,ylen,tlen*2));
%         c = cat(2,zeros((xlen*2),ylen/2,tlen*2),c);
%         c = cat(2,c,zeros((xlen*2),ylen/2,tlen*2));
        
        % determined the average rotation (time sample) rate using one
        % pixel location
        Tsampvec = datevec(mean(diff(cdnum + datenum(squeeze(t(xlen,ylen,:))/(24*60*60*1e6)))));
        Tsamp = Tsampvec(6);
        
        % set up inverted temporal axes, these are in radian frequency (rad/s)
        ws = 2*pi/Tsamp;
        wres = ws/(tlen*2);
        w = (wres*[-(tlen*2)/2:((tlen*2)/2)-1]);
        pw = w((tlen*2)/2+1:end);
        % low frequency cut off index
        wcutp = find(pw > wcutoff, 1 );
     
        % integrated polar PSD array
        ippsd = zeros(length(pw), length(kmag1));
        
        disp('FFT');
        c(isnan(c)) = 0;
        dspectrum = fftshift(fftn(c))./sqrt((xlen*2)*(ylen*2)*(tlen*2));
        % use positive frequency spectrum only (negative frequency spectrum
        % is complex conjugate of positive)
        pdspectrum = conj(flip(flip(flip(dspectrum(:,:,1:tlen),1),2),3));
        
        % compute power spectral density
        psd2 = 2*pdspectrum.*conj(pdspectrum);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Determine Doppler shifted dispersion relationship, due to current
        % for each positive frequency, change (kx,ky,w) coordinates to polar (|k|,theta,w)
        % integrate over theta to get a two dimensional power spectral density S(|k|,w)
        % This is the tricky part of unfolding the aliased spectrum
        disp('Doppler Shift');
   
        ppsd = zeros(size(xp));
        
        % PARFOR here, converts 3D PSD to polar coordinates, integrates out
        % theta to make a 2D spectrum
        for n=1:length(pw)
            ppsd = interp2(kxout,kyout,psd2(:,:,n),xp,yp,'linear',max(kx));
            ippsd(n,:) = sum(abs(ppsd),2)/size(ppsd,2)';
        end
        
        % unwrap the spectrum (it's aliased due to 10 RPM spin rate)
        % there are two unfolding points, so stack up the frequency vector
        % and the PSD
        pw2 = cat(2,pw,pw(end)+(wres.*(1:tlen)));
        pw3 = cat(2,pw2,pw2(end)+(wres.*(1:tlen)));
        ippsd2 = cat(1,ippsd,flipud(ippsd));
        ippsd3 = cat(1,ippsd2,ippsd);
        
        % plot aliased spectrum
        if ploton
            figure(3)
            contourf(kmag1,pw3,10*log10(abs(ippsd3)),30);
%             surf(kmag1,pw3,10*log10(abs(ippsd3)));
            axis([0 max(kx)/2 0 max(pw3)]);
            xlabel('|k| (rad/m)');
            ylabel('omega (rad/s)');
%             set(gca,'YTick',pw);
%             set(gca,'XTick',kmag1);
%             set(gca,'YTick',[0 pw(length(pw)/4) pw(length(pw)/2) pw(3*length(pw)/4) pw(end)]);
            colorbar;
            xlim([0 0.2]);
            caxis([0 50]);
            hold on;
            plot(kmag1,dispersionw,'-k');

            plot(kmag1,dispersionw - (kmag1.*cmax),'-k');
            plot(kmag1,dispersionw + (kmag1.*cmax),'-k');

            line([kmag1(kcutp1) kmag1(kcutp1)], [pw(1) pw(end)], 'Color', 'k');
            line([kmag1(kcutp2) kmag1(kcutp2)], [pw(1) pw(end)], 'Color', 'k');
            line([kmag1(1) kmag1(end)], [pw(wcutp) pw(wcutp)], 'Color', 'k');
        end
        
        % this reduces the search area on the spectrum within +/- the expected
        % maximum current parameter
        dmask = zeros(size(ippsd3));
        for iw = wcutp:length(pw3)
%                 dmask(iw,:) = (pw3(iw) >= dispersionw - (kmag1.*cmax) & ...
%                                (pw3(iw) <= dispersionw + (kmag1.*cmax)));
                dmask(iw,:) = (pw3(iw) >= dispersionw - 4*wres & ...
                               (pw3(iw) <= dispersionw + 4*wres));
        end
        dippsd = ippsd3.*dmask;
        dippsd(dippsd==0) = NaN;
        
        % erase part of  the aliased spectrum before filtering
        % this part needs to be adjusted depending on the spectrum
        
        % first aliasing point
        ka1 = pw3(tlen)^2./g;
        ka1idx = find(kmag1>=ka1,1);
        % mark the first aliasing point
        if ploton
            hold on
            line(kmag1,pw3(tlen)*ones(size(kmag1)),'color','r');
            line(kmag1,pw3(tlen+1)*ones(size(kmag1)),'color','r');
            line(kmag1,pw3(tlen*2)*ones(size(kmag1)),'color','r');
            line(kmag1,pw3(tlen*2+1)*ones(size(kmag1)),'color','r');
            hold off
        end
        % second aliasing point
        ka2 = pw3(2*tlen)^2./g;
        ka2idx = find(kmag1>=ka2,1);
        % mark the second aliasing point
        if ploton
            hold on
            plot(kmag1(ka2idx),pw3(tlen*2),'or')
            hold off
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % erase unusable spectrum
        % low frequency, high wavenumber, ENABLE/ADJUST THIS
        dippsd(tlen-1:tlen+1,:)=NaN;

        % mid-frequency band, all wavenumbers ENABLE/ADJUST THIS
        dippsd(2*tlen-1:2*tlen+1,:) = NaN;
%         
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
            colorbar;
            xlim([0 0.2]);
            caxis([0 20]);
            hold on;
            plot(kmag1,dispersionw,'-k');
            plot(kmag1,dispersionw - (kmag1.*cmax),'-k');
            plot(kmag1,dispersionw + (kmag1.*cmax),'-k');
            hold off
        end
       
        % plot red asterisks that will be used for the curve fit, make sure
        % they are on/near the shifted dispersion energy and not in the
        % noise or aliased spectrum, do it in three pieces
        %
        % 1) low cutoff to first aliasing point
        [wmax1,iwmax1] = max(dippsd(wcutp:tlen,kcutp1:ka1idx));
        iwmax1 = iwmax1+wcutp-1;  
        % 2) from first aliasing point to the second aliasing point
        [wmax2,iwmax2] = max(dippsd(tlen+2:tlen*2,ka1idx+1:ka2idx-1));
        iwmax2 = iwmax2+tlen+1;
        % 3) from second aliasing point to upper cutoff
        [wmax3,iwmax3] = max(dippsd(tlen*2+2:end,ka2idx+1:kcutp2));
        iwmax3 = iwmax3+tlen*2+1; 
        wmax = [wmax1 wmax2 wmax3];
        iwmax = [iwmax1 iwmax2 iwmax3];
        ik = [kcutp1:ka1idx ka1idx+1:ka2idx-1 ka2idx+1:kcutp2];
        if ploton
            hold on
            plot(kmag1(ik),pw3(iwmax),'r*');
            hold off
        end
        
        A=[dispersionw(ik).' kmag1(ik).'];
        b = pw3(iwmax).';
            
        inkmag = kmag1(ik);
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
        Ps = mean(wmax);
        [Pwmax, iPwmax] = max(wmax);
        % average noise region of S(|k|,w) to find average noise power spectral
        % density
        Pn = mean(mean(ippsd(wcutp:end,length(kmag1)/2:end)));
        if ploton
            figure(3)
            hold on
            line([kmag1(length(kmag1)/2), kmag1(length(kmag1)/2), kmag1(end), kmag1(end), kmag1(length(kmag1)/2)],...
                [pw(wcutp), pw(end), pw(end), pw(wcutp), pw(wcutp)],...
                'Color','k');
            hold off
        end

        % maximum 2D signal power is Pwmax at frequency pw3(iwmax(iPwmax)))
        iwmax = iwmax(iPwmax);
        if ploton
            figure(3)
            hold on
            plot(kmag1(ik(iPwmax)),pw3(iwmax),'or')
            title(['Maximum Signal Power ' num2str(10*log10(Pwmax)) ' dB at w = ' num2str(pw3(iwmax))]); 
            drawnow;
        end
        % now we know where maximum SNR is in the 2D spectrum, but need to
        % find it in the 3D spectrum as well, it will be at the same
        % frequency (w)
        
        % SET A BREAKPOINT HERE TO EXAMINE THE SPECTRUM AND ADJUST
        % PARAMETERS
    
        disp('Filter');
        % 3-D brick wall band pass filter
        % around the dispersion relationship, shifted for current,
        % +/-2 frequency bins, set all the rest = 0
        
        % make triple w spectrum
        temp = pdspectrum;
        pdspectrum = cat(3,pdspectrum,flip(flip(flip(temp,1),2),3));
        pdspectrum = cat(3,pdspectrum,temp);
        psd2 = pdspectrum.*conj(pdspectrum);
        
        % 3D signal power at previously found frequency
        figure(1)
        pcolor(kxout,kyout,10*log10(psd2(:,:,iwmax)));
        shading flat
        axis square
        caxis([0 50]);
        zoom(8)
        [Smax,ikxmax,ikymax] = max2d(psd2(:,:,iwmax));
        hold on
        plot(kxout(ikxmax,ikymax),kyout(ikxmax,ikymax),'or');
        title(['Maximum power at w=' num2str(pw3(iwmax),'%4.3f'),...
                               ' kx=' num2str(kxout(ikxmax,ikymax),'%4.3f'),...
                               ' ky=' num2str(kyout(ikxmax,ikymax),'%4.3f'),...
                               ' Dp=' num2str(mod(90-(atan2(kyout(ikxmax,ikymax),kxout(ikxmax,ikymax))*180/pi)+thetabox+dthetabox,360))])
        hold off
        % maximum 3D signal power is Smax, at (ikxmax,ikymax,iwmax)
        % set signal power threshold 
        threshold = Smax*0.1;
        
        Ax = reshape([reshape(dispersion2w,[],1) reshape(kmag2,[],1)]*x, (xlen*2), (ylen*2));
        % upper positive w limit, dispersion+current
        upwdc = Ax + 0.1*wres;
        % lower positive w limit, dispersion-current
        lpwdc = Ax - 0.1*wres;
        
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
                      kyout>0 &...
                      psd2(:,:,iw) >= threshold);
                      
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
            plot3(kxin,kyin,win,'.k');
            hold on
            mesh(KXFIT,KYFIT,WFIT);
            grid on
            hold off
            title(['Current Estimate: Ux = ' num2str(u(1),'%4.2f') '  Uy = ' num2str(u(2),'%4.2f')]);
            axis([-0.2 0.2 -0.2 0.2 0 max(pw3)])
        end 
       
        % now do it again, opening up filter width in frequency for waves
%         threshold = Smax*0.05;
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
                      kyout>0.001 &... 
                      psd2(:,:,iw) >= threshold);
%                           kxout>=-0.08 & kxout<=0.08 &...
%                           kyout>0.001 &...                    
            mask(:,:,iw) = gather(tmpmask);

            % use these points for wave inversion
            if sum(sum(mask(:,:,iw))) ~= 0
                kxin = [kxin; kxout(logical(mask(:,:,iw)))];
                kyin = [kyin; kyout(logical(mask(:,:,iw)))];
                win = [win; pw3(iw).*ones(sum(sum(mask(:,:,iw))),1)];
            end
        end
        % eliminate the 2nd alias frequencies (replicates DC noise)
        mask(:,:,tlen*2) = 0;
        mask(:,:,tlen*2+1) = 0;
        aidx = (win~=pw3(tlen*2)) & (win~=pw3(tlen*2+1));
        kxin = kxin(aidx);
        kyin = kyin(aidx);
        win = win(aidx);
        
        if ploton
            kxfit = min(kxin):0.001:max(kxin);
            kyfit = min(kyin):0.001:max(kyin);
            [KXFIT, KYFIT] = meshgrid(kxfit,kyfit);
            WFIT = sqrt(g*sqrt(KXFIT.^2+KYFIT.^2).*tanh(sqrt(KXFIT.^2+KYFIT.^2).*h)) + u(1).*KXFIT + u(2).*KYFIT;
            figure(6)
            plot3(kxin,kyin,win,'.k');
            axis([min2d(kxout) max2d(kxout) min2d(kyout) max2d(kyout) pw3(1) pw3(end)]);
            xlabel('kx');
            ylabel('ky');
            hold on
            mesh(KXFIT,KYFIT,WFIT);
            grid on
            hold off
        end
        
        pdspectrum = pdspectrum .* mask;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % plot filtered spectrum
        % for each positive frequency, change (kx,ky,w) coordinates to polar (|k|,theta,w)
        % integrate over theta to get a two dimensional power spectral density S(|k|,w)
        psd2 = pdspectrum.*conj(pdspectrum); 
        disp('Filtered Spectrum');
        
        % PARFOR, convert to polar, integrate out theta to make 2D spectrum
        for n=1:length(pw3)
            ppsd = interp2(kxout,kyout,psd2(:,:,n),xp,yp,'linear',max(kx));
            ippsd(n,:) = sum(abs(ppsd),2)/size(ppsd,2)';
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
        end
 
        % take the current out, don't currently use this...just checking to
        % see that (kxin,kyin,wun) is on the dispersion relationship, and
        % it is! 
        wun = win - [kxin kyin]*[u(1); u(2)];

        disp('IFFT');
        % create conjugate transpose for negative spectrum, concatenate
        ndspectrum = conj(flipdim(flipdim(flipdim(pdspectrum,1),2),3));
        ndspectrum = cat(3,zeros(xlen*2,ylen*2),ndspectrum);
        dspectrum = cat(3,ndspectrum,pdspectrum(:,:,1:end-1));
        waves = ifftn(ifftshift(dspectrum),'symmetric'); 
        zout = real(waves(xlen/2+1:xlen/2+xlen,ylen/2+1:ylen/2+ylen,tlen*2+tlen/2+1:tlen*2+tlen/2+tlen));
       
        % scale to Hs parameter, if you don't have it, just use Hs=1
        % uniform scaling
        sigma = std(zout,0,3);
%         zout1 = zout(:,:,tlen/2+1).*(Hs./(4*sigma));
        % relative scaling
        scale = sigma./max2d(sigma);
        zout = zout(:,:,tlen/2+1).*(Hs./(4*sigma)).*scale;
       
        sf=Hs./(4*std(reshape(zout,[],1)));
        zout = zout.*sf;

        tout = t(:,:,tlen/2+1);
        dnumout = cdnum(tlen/2+1);
        
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
%         end
       
        fdate = datestr(dnumout,'yyyymmddHHMMSS');

%         save([wavepath 'fs' fdate '.mat'],'kxout','kyout','w','dspectrum', 'A', 'b', 'Tsamp' );
        % save the (kx,ky,w) triplets for Doppler velocity inversion
        save([wavepath 'wkxky' fdate '.mat'], 'win','kxin','kyin','u');
        % save the waves
%         save([wavepath 'waves' fdate '.mat'],'zout','tout','dnumout');

    end
    toc
end

