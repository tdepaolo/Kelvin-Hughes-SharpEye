% read in Kelvin Hughes Sharp Eye radar data, generate single Frequency 1 and 
% Frequency 2 scans, ECEF (Earth centric Earth fixed) registered with ship data (lat, lon, heading)
% inputs: Pentek file name, scan time, file byte index (from file index)
% outputs: scan times, two frequency scans, fileheader, minimum and maximum turning
% indexes, squint, range vector, theta vector, polar lat/lon grids


function [scandnum,...
          f1time,f1mag,f1phi,...
          f2time,f2mag,f2phi,...
          fileheader,...
          minf1turning,maxf1turning,...
          minf2turning,maxf2turning,...
          squint,...
          r,theta,...
          lonp,latp] = read_SharpEye_sector(filename,byte_index)

    % get the parameters
    load('parameters.mat');
 
    % get file information
    load([datapath filename(16:21) '_info.mat']);
    
    % open the file               
    fid = fopen(filename);
    % fast forward to desired scan
    fseek(fid,byte_index,-1);

    % read Pentek Talon header and file header
%     talonheader = readTalonheader(fid);
    
    % read SharpEye file header
%     [fileheader, ferror] = readFileheader(fid);

    % convert to single, double
    sampletime = double(fileheader.sampletime);
    p1start = single(fileheader.p1start);
    p2start = double(fileheader.p2start);
    p3start = double(fileheader.p3start);           %Set to max if operating in two pulse mode
    fsamples = single(fileheader.fsamples);
    p1recst = single(fileheader.p1recst);
    p2recst = single(fileheader.p2recst);
    p3recst = single(fileheader.p3recst);
    recsamples = single(fileheader.recsamples);
    p1f1rng0 = single(fileheader.p1f1rng0);
    p2f1rng0 = single(fileheader.p2f1rng0);
    p3f1rng0 = single(fileheader.p3f1rng0);
    p1f2rng0 = single(fileheader.p1f2rng0);
    p2f2rng0 = single(fileheader.p2f2rng0);
    p3f2rng0 = single(fileheader.p3f2rng0);

    % sample numbers are 0-indexed (0-9999)
    % range bins are 1-indexed (MATLAB convention)
    % figure out range bins for each pulse
    % range bin 1 is 0 NM (at the radar)
    if (fileheader.p3start == MAX)
        p1samples = p2start - p1recst;          % number of samples
        p2samples = fsamples - p2recst;         % number of samples
        p3samples = NaN;
        p1f1ridx1 = p1recst - p1f1rng0 + 1;     % range array index, MATLAB convention
        p2f1ridx1 = p2recst - p2f1rng0 + 1;     
        p3f1ridx1 = NaN;
        p1f2ridx1 = p1recst - p1f2rng0;         % range index, KH convention
        p2f2ridx1 = p2recst - p2f2rng0;         % range index, KH convention
        p3f2ridx1 = NaN;
        % need these to index into sample array to get range 0
        p1f2radjust = -p1f2ridx1+1;             % range array index, MATLAB convention
        p2f2radjust = -p2f2ridx1+1;             % range array index, MATLAB convention
    else
        p1samples = p2start - p1recst;
        p2samples = p3start - p2recst;
        p3samples = fsamples - p3recst;
        p1f1ridx1 = p1recst - p1f1rng0 + 1;
        p2f1ridx1 = p2recst - p2f1rng0 + 1;
        p3f1ridx1 = p3recst - p3f1rng0 + 1;
        p1f2ridx1 = p1recst - p1f2rng0 + 1;
        p2f2ridx1 = p2recst - p2f2rng0 + 1;
        p3f2ridx1 = p3recst - p3f2rng0 + 1;
        % need to adjust these to range 0 sample
        p1f2radjust = -p1f2ridx1+1;
        p2f2radjust = -p2f2ridx1+1;
        p3f2radjust = -p3f2ridx1+1;
    end
    % maximum range sample recorded
    rsamples = p2samples + p2f1ridx1 - 1;

%     % check for lack of fileheader, then find the next frameheader
%     if ferror
%         nn = 0;
%         frameheader.framesync = fread(fid,1,'*uint32',0,'b');
%         while frameheader.framesync ~= FRAMESYNC
%             frameheader.framesync = fread(fid,1,'*uint32',0,'b');
%             dec2hex(frameheader.framesync,8);
%             nn = nn+1;
%         end
%         frameheader(1).framesync = frameheader.framesync;
%     else
        
%     end

    % read first frame header
    frameheader(1).framesync = fread(fid,1,'*uint32',0,'b');
    if frameheader(1).framesync ~= hex2dec('55AA0001');
        disp('framesync error');
        return
    end
    frameheader.frtimesec = fread(fid,1,'*uint32',0,'b');
    frameheader.frtimeusec = fread(fid,1,'*uint32',0,'b');
    frameheader.lastsmpls = fread(fid,1,'*uint32',0,'b');
    frameheader.turning = fread(fid,1,'*uint32',0,'b');
    frameheader.reserved1 = fread(fid,1,'*uint32',0,'b');
    frameheader.reserved2 = fread(fid,1,'*uint32',0,'b');
    frameheader.reserved3 = fread(fid,1,'*uint32',0,'b');

    % compute turning indexes, squint
    f1turning = bitand(frameheader(1).turning,F1TURNING,'uint32');
    f2turning = bitshift(bitand(frameheader(1).turning,F2TURNING,'uint32'),-16);
    squint = mod(single(f2turning) - single(f1turning),4096);
    % get the time stamps for this scan
    frtimesec = double(frameheader(1).frtimesec);
    scandnum = SharpEyeTime2UTC(frtimesec);

    % there are twice as many samples since there are two frequencies, there
    % may be less samples since frames are variable length
    sample = zeros(1,fileheader.recsamples*2+frameheaderlen,'uint32');

    % maximum number of pulses per turning (determined empirically by looking
    % at the recordings)
    if isnan(p3samples)
        pulses_per_turning = pulses_per_turning_12NM;
        % maximum number of pulses per range bin, due to pulses overlapping in range
        pulses_per_range_mode = pulses_per_range_mode_12NM;
    else
        pulses_per_turning = pulses_per_turning_24NM;
        pulses_per_range_mode = pulses_per_range_mode_24NM;
    end

    % frame and pulse counter per turning index
    f1pulsecount = 1;
    f2pulsecount = 1;
    extrasamples = 0;

    %load in the ShipData file for the date of interest
%     [sdnum, slat, slon, scog, ssog, shead] = load_shipData(datestr(scandnum,'yyyymmdd'));

    f1saveon = 1;
    f2saveon = 1;
    F2DONE = 0;
    F1DONE = 0;
    
    % set up scan array
    % indexing of the scan arrays is (range,azimuth,pulses)
    % Pulses are stored in the order received for the frame
    % in 6NM mode it's P1/P2/P1/P2... in the pulses dimension

    %                 range vector in radar coordinates
    r = ((0:rsamples-1)'.*rres);
    %                 azimuth vector in compass coordinates
    % need to check for spanning the heading line (turning index 0)
    if (maxf2turning>minf2turning)
        theta = mod((single(minf1turning):single(maxf2turning)).*(360./4096) + single(HL_offset),360).*pi/180;
    else
        sector_width = 4095 - minf1turning + 1 + maxf2turning;
        theta =  mod((single(minf1turning):single(minf1turning)+sector_width).*(360./4096) + single(HL_offset),360).*pi/180;
    end
    
    f1mag =  NaN.*ones(rsamples,length(theta),pulses_per_turning*pulses_per_range_mode,'single');
    f1phi =  NaN.*ones(rsamples,length(theta),pulses_per_turning*pulses_per_range_mode,'single');
    f1time = NaN.*ones(rsamples,length(theta),pulses_per_turning*pulses_per_range_mode,'double');
    f2mag =  NaN.*ones(rsamples,length(theta),pulses_per_turning*pulses_per_range_mode,'single');
    f2phi =  NaN.*ones(rsamples,length(theta),pulses_per_turning*pulses_per_range_mode,'single');
    f2time = NaN.*ones(rsamples,length(theta),pulses_per_turning*pulses_per_range_mode,'double');

    % find range and azimuth indexes to Doppler process, defines a wedge, box is inside of the wedge, from parameters
    r1 = find(r>=rmin,1,'first');
    r2 = find(r>=rmax,1,'first');
    [val, az1] = min(abs(theta-azmin*pi/180));
    [val, az2] = min(abs(theta-azmax*pi/180));

    %  azimuth array in Cartesian coordinates for MATLAB sin,cos
    mtheta = pi/2 - theta;
    % radar scan grid, e,n from radar origin
    e = r*cos(mtheta);
    n = r*sin(mtheta);
    % registered polar grid
    [lonp, latp] = km2lonlat(khlongitude, khlatitude, e./1000, n./1000);

    st = sprintf('read_SharpEye: minf1turning %d',minf1turning);
    disp(st);
    st = sprintf('read_SharpEye: maxf1turning %d',maxf1turning);
    disp(st);
    st = sprintf('read_SharpEye: minf2turning %d',minf2turning);
    disp(st);        
    st = sprintf('read_SharpEye: maxf2turning %d',maxf2turning);
    disp(st);

    while ~F1DONE

        % compute turning indexes
        f1turning = bitand(frameheader(1).turning,F1TURNING,'uint32');
        f2turning = bitshift(bitand(frameheader(1).turning,F2TURNING,'uint32'),-16);  
        
        % frames are variable size, and you don't know how many samples were
        % actually written in the last frame until you read the next frame
        % header. At this point, the frame header for the next frame has been
        % read (except for the first frame of the file), and is at the end of
        % the "sample" array, and any additional samples for the next frame.

        % From the last sample array, grab the extra samples from the next frame if 
        % there are any, and put them at the start of the next sample array
        % (not applicable in the first frame of the recording)
        if extrasamples ~=0
            sample(1:extrasamples) = sample(end-extrasamples+1:end);
        end

        % Now read all data samples in this frame, and the
        % next frame header, and perhaps some extra samples from the next frame

        sample(1+extrasamples:end) = fread(fid,fileheader.recsamples*2+frameheaderlen-extrasamples,'*uint32',0,'b');
        findex = find(sample==FRAMESYNC,1,'first');
        frameheader(2).framesync = sample(findex);
        frameheader(2).frtimesec = sample(findex+1);
        frameheader(2).frtimeusec = sample(findex+2);
        frameheader(2).lastsmpls = sample(findex+3);
        frameheader(2).turning = sample(findex+4);
        frameheader(2).reserved1 = sample(findex+5);
        frameheader(2).reserved2 = sample(findex+6);
        frameheader(2).reserved3 = sample(findex+7);

        extrasamples = fileheader.recsamples*2 - frameheader(2).lastsmpls;

        % compute turning angles for frame 2
        nextf1turning = bitand(frameheader(2).turning,F1TURNING,'uint32');
        nextf2turning = bitshift(bitand(frameheader(2).turning,F2TURNING,'uint32'),-16);


        % FREQUENCY 1
        % magnitude and phase samples for frame 1, fill end with zeros if a short frame
        phase = [uint16(bitand(sample(1:2:frameheader(2).lastsmpls-1),PHASE,'uint32')) zeros(1,fileheader.recsamples-frameheader(2).lastsmpls/2,'uint16')];
        signbit = bitget(phase,15,'uint16');
        % assume phase is two's complement, make it 16 bits
        phase(signbit==1) = bitor(phase(signbit==1),SIGNBIT,'uint16');          % extend the sign bit
        phase(signbit==1) = bitcmp(phase(signbit==1))+1;                        % complement, add 1
        phase = int16(phase);                                                   % unsigned to signed
        phase(signbit==1) = -phase(signbit==1);                                 % negate
        phase = single(phase);                                                  % float it
        % mask of the magnitude, shift it down, convert to float
        magnitude = [single(bitshift(bitand(sample(1:2:frameheader(2).lastsmpls-1),MAGNITUDE,'uint32'),-15)) zeros(1,fileheader.recsamples-frameheader(2).lastsmpls/2,'single')];


        % translate sample number to range bin, frequency 1 samples start after range 0
        % pulse 1
        f1mag(p1f1ridx1:p1f1ridx1+p1samples-1,mod(single(f1turning)-single(minf1turning)+1,4096),f1pulsecount) = ...
            magnitude(1:p1samples);
        f1phi(p1f1ridx1:p1f1ridx1+p1samples-1,mod(single(f1turning)-single(minf1turning)+1,4096),f1pulsecount) = ...
            phase(1:p1samples);
        % compute timestamps in usec
        f1time(p1f1ridx1:p1f1ridx1+p1samples-1,mod(single(f1turning)-single(minf1turning)+1,4096),f1pulsecount) = ... 
                (double(frameheader(1).frtimesec) - frtimesec).*1.0e6 + ...
                 double(frameheader(1).frtimeusec) + ...
                (double(p1f1ridx1:p1f1ridx1+p1samples-1)-1).*sampletime.*1.0e-3;

        f1pulsecount = f1pulsecount + 1;


        % pulse 2
        f1mag(p2f1ridx1:p2f1ridx1+p2samples-1,mod(single(f1turning)-single(minf1turning)+1,4096),f1pulsecount) = ...
            magnitude(p1samples+1:p1samples+p2samples);
        f1phi(p2f1ridx1:p2f1ridx1+p2samples-1,mod(single(f1turning)-single(minf1turning)+1,4096),f1pulsecount) = ...
            phase(p1samples+1:p1samples+p2samples);
        % compute timestamps in usec
        f1time(p2f1ridx1:p2f1ridx1+p2samples-1,mod(single(f1turning)-single(minf1turning)+1,4096),f1pulsecount) = ... 
                (double(frameheader(1).frtimesec) - frtimesec).*1.0e6 + ...
                 double(frameheader(1).frtimeusec) + (p2start.*sampletime.*1.0e-3) + ...
                (double(p2f1ridx1:p2f1ridx1+p2samples-1)-1).*sampletime.*1.0e-3;

        f1pulsecount = f1pulsecount + 1;

        % pulse 3
        if ~isnan(p3samples)
            f1mag(p3f1ridx1:p3f1ridx1+p3samples-1,mod(single(f1turning)-single(minf1turning)+1,4096),f1pulsecount) = ...
                magnitude(p1samples+p2samples+1:p1samples+p2samples+p3samples);
            f1phi(p3f1ridx1:p3f1ridx1+p3samples-1,mod(single(f1turning)-single(minf1turning)+1,4096),f1pulsecount) = ...
                phase(p1samples+p2samples+1:p1samples+p2samples+p3samples);
            % compute timestamps in usec
            f1time(p3f1ridx1:p3f1ridx1+p3samples-1,mod(single(f1turning)-single(minf1turning)+1,4096),f1pulsecount) = ... 
                (double(frameheader(1).frtimesec) - frtimesec).*1.0e6 + ...
                 double(frameheader(1).frtimeusec) + (p3start.*sampletime.*1.0e-3) + ...
                (double(p3f1ridx1:p3f1ridx1+p3samples-1)-1).*sampletime.*1.0e-3;

            f1pulsecount = f1pulsecount + 1;
        end

        if f1turning ~= nextf1turning
            % check for turning index adjustment
%                 if (single(f1turning)-single(nextf1turning))>1
%                     minf1turning = nextf1turning;
%                     maxf1turning = f1turning;
%                 end

            if f1turning == maxf1turning
                disp('read_SharpEye: end of F1 scan');
                if f1saveon
                    F1DONE = 1;
                else
                    disp('read_SharpEye: skipping end of first F1 scan');
                    % reinitialize arrays
                    f1mag(:,:,:) = NaN;
                    f1phi(:,:,:) = NaN;
                    f1time(:,:,:) = NaN;
                    f1saveon = 1;
                end
            end
            f1pulsecount = 1;
        end

        % FREQUENCY 2
        if ~F2DONE
            phase = [uint16(bitand(sample(2:2:frameheader(2).lastsmpls),PHASE,'uint32')) zeros(1,fileheader.recsamples-frameheader(2).lastsmpls/2,'uint16')];
            signbit = bitget(phase,15,'uint16');
            % assume phase is two's complement, make it 16 bits
            phase(signbit==1) = bitor(phase(signbit==1),SIGNBIT,'uint16');          % extend the sign bit
            phase(signbit==1) = bitcmp(phase(signbit==1))+1;                        % complement, add 1
            phase = int16(phase);                                                   % unsigned to signed
            phase(signbit==1) = -phase(signbit==1);                                 % negate
            phase = single(phase);                                                  % float it
            % mask of the magnitude, shift it down, convert to float
            magnitude = [single(bitshift(bitand(sample(2:2:frameheader(2).lastsmpls),MAGNITUDE,'uint32'),-15)) zeros(1,fileheader.recsamples-frameheader(2).lastsmpls/2,'single')];

            % old way, no range correction for F2
            f2mag(1:p1f2ridx1+p1samples,mod(single(f2turning)-single(minf1turning)+1,4096),f2pulsecount) = ...
                magnitude(p1f2radjust:p1samples);
            f2phi(1:p1f2ridx1+p1samples,mod(single(f2turning)-single(minf1turning)+1,4096),f2pulsecount) = ...
                phase(p1f2radjust:p1samples);
            % compute timestamps in usec
            f2time(1:p1f2ridx1+p1samples,mod(single(f2turning)-single(minf1turning)+1,4096),f2pulsecount) = ... 
                    (double(frameheader(1).frtimesec) - frtimesec).*1.0e6 + ...
                     double(frameheader(1).frtimeusec) + ...
                    (double(1:p1f2ridx1+p1samples)-1).*sampletime.*1.0e-3;

            f2pulsecount = f2pulsecount + 1;

            % pulse 2
            f2mag(1:p2f2ridx1+p2samples,mod(single(f2turning)-single(minf1turning)+1,4096),f2pulsecount) = ...
                magnitude(p1samples+p2f2radjust:p1samples+p2samples);
            f2phi(1:p2f2ridx1+p2samples,mod(single(f2turning)-single(minf1turning)+1,4096),f2pulsecount) = ...
                phase(p1samples+p2f2radjust:p1samples+p2samples);
            % compute timestamps in usec
            f2time(1:p2f2ridx1+p2samples,mod(single(f2turning)-single(minf1turning)+1,4096),f2pulsecount) = ... 
                    (double(frameheader(1).frtimesec) - frtimesec).*1.0e6 + ...
                     double(frameheader(1).frtimeusec) + (p2start.*sampletime.*1.0e-3) + ...
                    (double(1:p2f2ridx1+p2samples)-1).*sampletime.*1.0e-3;

            f2pulsecount = f2pulsecount + 1;

            % pulse 3
            if ~isnan(p3samples)
                f2mag(1:p3f2ridx1+p3samples,mod(single(f2turning)-single(minf1turning)+1,4096),f2pulsecount) = ...
                    magnitude(p1samples+p2samples+p3f2radjust:p1samples+p2samples+p3samples);
                f2phi(1:p3f2ridx1+p3samples,mod(single(f2turning)-single(minf1turning)+1,4096),f2pulsecount) = ...
                    phase(p1samples+p2samples+p3f2radjust:p1samples+p2samples+p3samples);
                % compute timestamps in usec
                f2time(1:p3f2ridx1+p3samples,mod(single(f2turning)-single(minf1turning)+1,4096),f2pulsecount) = ... 
                    (double(frameheader(1).frtimesec) - frtimesec).*1.0e6 + ...
                     double(frameheader(1).frtimeusec) + (p3start.*sampletime.*1.0e-3) + ...
                    (double(1:p3f2ridx1+p3samples)-1).*sampletime.*1.0e-3;

                f2pulsecount = f2pulsecount + 1;
            end

%             % translate sample number to range bin, frequency 2 samples start
%             % BEFORE range 0 (garbage), so those samples are thrown away
%             % there is a 2 range bin offset for F2, don't know
%             % why....adjusting it here.
%             % pulse 1
%             f2mag(1:p1f2ridx1+p1samples-2,mod(single(f2turning)-single(minf1turning)+1,4096),f2pulsecount) = ...
%                 magnitude(p1f2radjust+2:p1samples);
%             f2phi(1:p1f2ridx1+p1samples-2,mod(single(f2turning)-single(minf1turning)+1,4096),f2pulsecount) = ...
%                 phase(p1f2radjust+2:p1samples);
%             % compute timestamps in usec
%             f2time(1:p1f2ridx1+p1samples,mod(single(f2turning)-single(minf1turning)+1,4096),f2pulsecount) = ... 
%                     (double(frameheader(1).frtimesec) - frtimesec).*1.0e6 + ...
%                      double(frameheader(1).frtimeusec) + ...
%                     (double(1:p1f2ridx1+p1samples)-1).*sampletime.*1.0e-3;
% 
%             f2pulsecount = f2pulsecount + 1;
% 
%             % pulse 2
%             f2mag(1:p2f2ridx1+p2samples-2,mod(single(f2turning)-single(minf1turning)+1,4096),f2pulsecount) = ...
%                 magnitude(p1samples+p2f2radjust+2:p1samples+p2samples);
%             f2phi(1:p2f2ridx1+p2samples-2,mod(single(f2turning)-single(minf1turning)+1,4096),f2pulsecount) = ...
%                 phase(p1samples+p2f2radjust+2:p1samples+p2samples);
%             % compute timestamps in usec
%             f2time(1:p2f2ridx1+p2samples,mod(single(f2turning)-single(minf1turning)+1,4096),f2pulsecount) = ... 
%                     (double(frameheader(1).frtimesec) - frtimesec).*1.0e6 + ...
%                      double(frameheader(1).frtimeusec) + (p2start.*sampletime.*1.0e-3) + ...
%                     (double(1:p2f2ridx1+p2samples)-1).*sampletime.*1.0e-3;
% 
%             f2pulsecount = f2pulsecount + 1;
% 
%             % pulse 3
%             if ~isnan(p3samples)
%                 f2mag(1:p3f2ridx1+p3samples-1,mod(single(f2turning)-single(minf1turning)+1,4096),f2pulsecount) = ...
%                     magnitude(p1samples+p2samples+p3f2radjust+1:p1samples+p2samples+p3samples);
%                 f2phi(1:p3f2ridx1+p3samples-1,mod(single(f2turning)-single(minf1turning)+1,4096),f2pulsecount) = ...
%                     phase(p1samples+p2samples+p3f2radjust+1:p1samples+p2samples+p3samples);
%                 % compute timestamps in usec
%                 f2time(1:p3f2ridx1+p3samples,mod(single(f2turning)-single(minf1turning)+1,4096),f2pulsecount) = ... 
%                     (double(frameheader(1).frtimesec) - frtimesec).*1.0e6 + ...
%                      double(frameheader(1).frtimeusec) + (p3start.*sampletime.*1.0e-3) + ...
%                     (double(1:p3f2ridx1+p3samples)-1).*sampletime.*1.0e-3;
% 
%                 f2pulsecount = f2pulsecount + 1;
%             end

            if f2turning ~= nextf2turning
                % check for turning index adjustment
%                 if (single(f2turning)-single(nextf2turning))>1
%                     minf2turning = nextf2turning;
%                     maxf2turning = f2turning;
%                 end
%                         figure(10)
%                         hold on
%                         t2 = reshape(f2time(:,f2turning+1,:),[],1);
%                         plot(t1,'g');

                % check for end of scan
                if f2turning == maxf2turning
                    disp('read_SharpEye: end of F2 scan');
                    if f2saveon
                        F2DONE = 1;
                    end
                end
                f2pulsecount = 1;
            end
        end
        % save frame 2 header as frame 1 for next pass
        frameheader(1) = frameheader(2);
    end

    fclose(fid);
end
              