function [Doppler_frequency_fft1, Doppler_frequency_pp1, conf1, Doppler_frequency_fft2, Doppler_frequency_pp2, conf2] = compute_Doppler_frequency_separate(f1mag,f1phi,f1t,f2mag,f2phi,f2t,dt)
    % inputs: magnitude, phase, time vectors, and sample time
    % outputs: Doppler frequency using FFT, Doppler frequency using pulse pair (pp)
  
    % turn plotting off/on
    ploton = 0;
    
    % set up frequency vectors
    fftlen = 64;
    % sampling rate
    fs = 1/(dt.*1e-6);
    % frequency resolution
    fres = fs/fftlen;
    % frequency vector
    f = (-fftlen/2:(fftlen/2)-1)*fres;
    % hanning window for fft
    hwin = hann(fftlen,'symmetric');

    % convert input to the fftlen for windowing
    if length(f1t)>fftlen
        f1t = f1t(1:fftlen);
        f1mag = f1mag(1:fftlen);
        f1phi = f1phi(1:fftlen);
    elseif length(f1t)<fftlen
        for tl=length(f1t)+1:fftlen
            f1t = cat(1,f1t,f1t(end)+dt);
        end
        f1mag = cat(1,f1mag,zeros(fftlen-length(f1mag),1));
        f1phi = cat(1,f1phi,zeros(fftlen-length(f1phi),1));
    end
    
    if length(f2t)>fftlen
        f2t = f2t(1:fftlen);
        f2mag = f2mag(1:fftlen);
        f2phi = f2phi(1:fftlen);
    elseif length(f2t)<fftlen
        for tl=length(f2t)+1:fftlen
            f2t = cat(1,f2t,f2t(end)+dt);
        end
        f2mag = cat(1,f2mag,zeros(fftlen-length(f2mag),1));
        f2phi = cat(1,f2phi,zeros(fftlen-length(f2phi),1));
    end
        
    % F1
    % unwrap the phase
    uphi = unwrap(f1phi);
        
    % make I,Q
    I = f1mag.*cos(uphi);
    Q = f1mag.*sin(uphi);

    % low pass filter
    navg = 4;
    b = ones(navg,1).*(1/navg);
    a = 1;

    smoothp = filtfilt(b,a,double(uphi));
    smoothI = filtfilt(b,a,double(I));
    smoothQ = filtfilt(b,a,double(Q));

    % Doppler spectrum   
    S = fftshift(fft(hwin.*(smoothI+1i.*smoothQ),fftlen));

    % generate power spectrum, interpolate to 1 Hz resolution
    ff = min(f):floor(max(f));
    PSD = (S.*conj(S))./fftlen;
    SS = interp1(f,PSD,ff);
    
    % compute first moment to get Doppler frequency
    Doppler_frequency_fft1 = sum(SS.*ff)./sum(SS);
    
    % compute SNR
%     N0 = mean(SS(3*fftlen/4:end));
%     snr1 = (sum(SS)-N0)./N0;

    % Zrnic pulse pair
    % take autocorrelation
    R = (1./length(smoothp))*sum(xcorr(conj(I+1i*Q),I+1i*Q));
    % pulse pair, determine slope (frequency) between each pair
    dpp = diff(smoothp./(2*pi))./(diff(f1t).*1e-6);
    conf1 = abs(sum(exp(1i*diff(smoothp))))./(length(smoothp)-1);
    % generate distribution (pdf)p2tin
    dpp_dist = histc(dpp,f);
%     dpp_dist = histcounts(dpp,[f fs/2],'Normalization','pdf');
    % spline fit,  distribution
    DPP_DIST = interp1(f,dpp_dist,ff);
    % compute first moment to get Doppler frequency
    Doppler_frequency_pp1 = sum(DPP_DIST.*ff)./sum(DPP_DIST);
    % take maximum
%     [val, idx] =  max(DPP_DIST);
%     Doppler_frequency_pp1 = ff(idx);

    if ploton
        figure(5)
        subplot(2,1,1)
        plot(f1t,I,'b',f1t,Q,'g',f1t,smoothI,'ob',f1t,smoothQ,'og');
%         ylim([-10 10]);
        subplot(2,1,2)
        semilogy(f,PSD,'b',ff,SS,'g');
%         xlim([-500 500]);
        title(['F1 Doppler Frequency = ' num2str(Doppler_frequency_fft1,'%4.2f')]);
        
        figure(6)
        subplot(2,1,1)
        plot(f1t,smoothp);
        subplot(2,1,2)
        plot(f,dpp_dist,'b',ff,DPP_DIST,'g');
        title(['F1 PP Frequency = ' num2str(Doppler_frequency_pp1,'%4.2f')]);
        drawnow;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % F2   
    % unwrap the phase
    uphi = unwrap(f2phi);
       
    % make I,Q
    I = f2mag.*cos(uphi);
    Q = f2mag.*sin(uphi);

    % low pass filter
    navg = 4;
    b = ones(navg,1).*(1/navg);
    a = 1;
    smoothp = filtfilt(b,a,double(uphi));
    smoothI = filtfilt(b,a,double(I));
    smoothQ = filtfilt(b,a,double(Q));

    % Doppler spectrum
    S = fftshift(fft(hwin.*(smoothI+1i.*smoothQ),fftlen));

    % generate power spectrum, interpolate to 1 Hz resolution
    ff = min(f):floor(max(f));
    PSD = (S.*conj(S))./fftlen;
    SS = interp1(f,PSD,ff);
    
    % compute first moment to get Doppler frequency
    Doppler_frequency_fft2 = sum(SS.*ff)./sum(SS);
    
    % compute SNR
%     N0 = mean(SS(3*fftlen/4:end));
%     snr2 = (sum(SS)-N0)./N0;
     % pulse pair, determine slope (frequency) between each pair
    dpp = diff(smoothp./(2*pi))./(diff(f2t).*1e-6);
%     conf2 = abs(sum(exp(1i*diff(smoothp))))./(length(smoothp)-1);
    conf2 = abs(sum(exp(1i*diff(smoothp))))./(length(smoothp)-1);
    % generate distribution (pdf)
%     dpp_dist = histcounts(dpp,[f fs/2],'Normalization','pdf');
    dpp_dist = histc(dpp,f);
    % spline fit,  distribution
    DPP_DIST = interp1(f,dpp_dist,ff);
    % compute first moment to get Doppler frequency
    Doppler_frequency_pp2 = sum(DPP_DIST.*ff)./sum(DPP_DIST);
    % take maximum
%     [val, idx] =  max(DPP_DIST);
%     Doppler_frequency_pp2 = ff(idx);
    
    if ploton
        figure(5)
        subplot(2,1,1)
        hold on
        plot(f2t,I,'r',f2t,Q,'m',f2t,smoothI,'or',f2t,smoothQ,'om');
%         ylim([-10 10])
        legend('F1I','F1Q','F2I','F2Q');
        hold off
        subplot(2,1,2)
        hold on
        semilogy(f,PSD,'r',ff,SS,'m');
        legend('',['F1 Doppler f=' num2str(Doppler_frequency_fft1)]','',['F2 Doppler f=' num2str(Doppler_frequency_fft2)]');
        drawnow;
        hold off
        
        figure(6)
        subplot(2,1,1)
        hold on
        plot(f1t,smoothp,'r');
        legend('F1 PP','F2 PP');
        hold off
        subplot(2,1,2)
        hold on
        plot(f,dpp_dist,'r',ff,DPP_DIST,'m');
        legend('',['F1 PP f=' num2str(Doppler_frequency_pp1)]','',['F2 PP f=' num2str(Doppler_frequency_pp2)]');
        drawnow;
        hold off
    end
    
%     % pulse pair, determine slope (frequency) between each pair
%     dpp = diff(smoothp./(2*pi))./(diff(f2t).*1e-6);
%     % generate distribution (pdf)
%     dpp_dist = histc(dpp,f);
%     % spline fit,  distribution
%     DPP_DIST = interp1(f,dpp_dist,ff);
%     % compute first moment to get Doppler frequency
%     Doppler_frequency_pp = sum(DPP_DIST.*ff)./sum(DPP_DIST);
%     % take maximum
%     [val, idx] =  max(DPP_DIST);
%     Doppler_frequency_pp = ff(idx);
    

    
end