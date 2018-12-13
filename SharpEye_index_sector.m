% index a Kelvin Hughes Sharp Eye Pentek file, generate an index file
% (ASCII) and a .mat file for use in reading in any scan in the file.
% This script does not generate and scans, frames, etc....only an index.
% The index marks the end of each F2 scan, which leads the F1 scans by
% squint turnings.  The next byte in the file after the index value is the
% start of the next F2 scan.

% Saves info on the scan:  fileheader, min/max turning inedexes

clear all
close all

% get the parameters
load('parameters.mat');

% pop up window to pick file to index
filename = uigetfile([datapath '*.dat']);

% open a file               
fid = fopen([datapath filename]);

% % read Pentek Talon header
talonheader=readTalonheader(fid);
% 
% % read SharpEye file header
[fileheader, ferror] = readFileheader(fid);
% 
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

% figure out range bins for each pulse
% range bin 1 is 0 NM (at the radar)
if (p3start == MAX)
    p1samples = p2start - p1recst;
    p2samples = recsamples + 1 - p2recst;
    p3samples = NaN;
    p1f1ridx1 = p1recst - p1f1rng0 + 1;
    p2f1ridx1 = p2recst - p2f1rng0 + 1;
    p3f1ridx1 = NaN;
    p1f2ridx1 = p1recst - p1f2rng0;
    p2f2ridx1 = p2recst - p2f2rng0;
    p3f2ridx1 = NaN;
    % need to adjust these to range 0 sample
    p1f2radjust = -p1f2ridx1+1;
    p2f2radjust = -p2f2ridx1+1;
else
    p1samples = p2start - p1recst;
    p2samples = p3start - p2recst;
    p3samples = recsamples + 1 - p3recst;
    p1f1ridx1 = p1recst - p1f1rng0 + 1;
    p2f1ridx1 = p2recst - p2f1rng0 + 1;
    p3f1ridx1 = p3recst - p3f1rng0 + 1;
    p1f2ridx1 = p1recst - p1f2rng0;
    p2f2ridx1 = p2recst - p2f2rng0;
    p3f2ridx1 = p3recst - p3f2rng0;
    % need to adjust these to range 0 sample
    p1f2radjust = -p1f2ridx1+1;
    p2f2radjust = -p2f2ridx1+1;
    p3f2radjust = -p3f2ridx1+1;
end
% maximum range sample recorded
rsamples = p2samples + p2f1ridx1 - 1;

% open index file
% index_fid = fopen([d(1).name(1:12) '_index.txt'],'w','n','US-ASCII');
% start an index matrix
index = [];

if ferror
    % find the frame header
    frameheader.framesync = fread(fid,1,'*uint32',0,'b');
    while frameheader.framesync ~= FRAMESYNC
        frameheader.framesync = fread(fid,1,'*uint32',0,'b');
%         dec2hex(frameheader.framesync,8)
    end
    frameheader(1).framesync = frameheader.framesync;
else
    % read first frame header
    frameheader(1).framesync = fread(fid,1,'*uint32',0,'b');
end

frameheader.frtimesec = fread(fid,1,'*uint32',0,'b');
frameheader.frtimeusec = fread(fid,1,'*uint32',0,'b');
frameheader.lastsmpls = fread(fid,1,'*uint32',0,'b');
frameheader.turning = fread(fid,1,'*uint32',0,'b');
frameheader.reserved1 = fread(fid,1,'*uint32',0,'b');
frameheader.reserved2 = fread(fid,1,'*uint32',0,'b');
frameheader.reserved3 = fread(fid,1,'*uint32',0,'b');

% compute turning indexes, squint
disp('First turning indexes in file');
f1turning = bitand(frameheader(1).turning,F1TURNING,'uint32');
st = sprintf('f1turning %d',f1turning);
disp(st);
f2turning = bitshift(bitand(frameheader(1).turning,F2TURNING,'uint32'),-16);
st = sprintf('f2turning %d',f2turning);
disp(st);
squint = f2turning - f1turning;

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
extrasamples = 0;

f1frtimesec = double(frameheader(1).frtimesec);
f1filednum = SharpEyeTime2UTC(f1frtimesec);
f2frtimesec = double(frameheader(1).frtimesec);
f2filednum = SharpEyeTime2UTC(f2frtimesec);

firstf1turning = f1turning;
framecount = 0;
f1_turning_indexes = zeros(4096,1);
f2_turning_indexes = zeros(4096,1);
SETUP = 0;
SEEK = 1;
SAVE = 2;
state = SETUP;
scancount = 0;

while ~feof(fid)
    try
        % compute turning indexes
        f1turning = bitand(frameheader(1).turning,F1TURNING,'uint32');
        f2turning = bitshift(bitand(frameheader(1).turning,F2TURNING,'uint32'),-16);

        % frames are variable size, and you don't know how many samples were
        % actually written in the last frame until you read the next frame
        % header.  From the second frame on, there may be extra samples at the
        % end of the array that may be from the next frame.

        % grab the extra samples from the last frame if there are any
        % (not applicable in the first frame of the recording)
        if extrasamples ~=0
            sample(1:extrasamples) = sample(end-extrasamples+1:end);
        end

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

        switch state
            case SETUP 
                % determine all turning indexes in a scan
                framecount = framecount + 1;
                f1_turning_indexes(f1turning+1) = 1;
                f2_turning_indexes(f2turning+1) = 1;
                if nextf1turning==firstf1turning && framecount>pulses_per_turning
                    if sum(f1_turning_indexes) == 4096     % full scan 
                        minf1turning = 0;
                        maxf1turning = 4095;
                        minf2turning = 0;
                        maxf2turning = 4095;
                    else                                % not full scan
                    % this is tricky, zero indexing vs. 1 indexing
                        zindex = find(f1_turning_indexes==0);
                        minf1turning = mod(zindex(end)+1,4096);
                        minf1turning = mod(minf1turning-1,4095);
                        maxf1turning = mod(zindex(1)-1,4096);
                        maxf1turning = mod(maxf1turning-1,4095);
                        zindex = find(f2_turning_indexes==0);
                        minf2turning = mod(zindex(end)+1,4096);
                        minf2turning = mod(minf2turning-1,4095);
                        maxf2turning = mod(zindex(1)-1,4096);
                        maxf2turning = mod(maxf2turning-1,4095);
%                         minf1turning = find(f1_turning_indexes,1,'first')-1;
%                         maxf1turning = find(f1_turning_indexes,1,'last')-1;
%                         minf2turning = find(f2_turning_indexes,1,'first')-1;
%                         maxf2turning = find(f2_turning_indexes,1,'last')-1;
                    end

                    f1frtimesec = double(frameheader(2).frtimesec);
                    f1filednum = SharpEyeTime2UTC(f1frtimesec);

                    % change state
                    st = sprintf('minf1turning %d',minf1turning);
                    disp(st);
                    st = sprintf('maxf1turning %d',maxf1turning);
                    disp(st);
                    st = sprintf('minf2turning %d',minf2turning);
                    disp(st);
                    st = sprintf('maxf2turning %d',maxf2turning);
                    disp(st);
                    state = SEEK;
                end
                    % end of SETUP

            case SEEK
                % spin here until we get to the start of the next scan

                if nextf2turning==minf2turning 
                    disp('starting F2 scan');
                    f1frtimesec = double(frameheader(1).frtimesec);
                    f1filednum = SharpEyeTime2UTC(f1frtimesec);
                    f2frtimesec = double(frameheader(1).frtimesec);
                    f2filednum = SharpEyeTime2UTC(f2frtimesec);
                    state = SAVE;
                end

            case SAVE
                if f1turning ~= nextf1turning
                    % check for end of scan
                    if f1turning == maxf1turning
                        f2frtimesec = double(frameheader(2).frtimesec);
                        f2filednum = SharpEyeTime2UTC(f2frtimesec);
                        st = sprintf('end of F1 scan %s',datestr(f2filednum,'yyyymmddHHMMSS'));
                        disp(st);
                        scancount = scancount + 1;
                    end
                end
                if f2turning ~= nextf2turning
                    % check for end of scan
                    if f2turning == maxf2turning
                        f2frtimesec = double(frameheader(2).frtimesec);
                        f2filednum = SharpEyeTime2UTC(f2frtimesec);
                        st = sprintf('end of F2 scan %s',datestr(f2filednum,'yyyymmddHHMMSS'));
                        disp(st);
%                         fprintf(index_fid,'date: %s - fileposition: %d\n',datestr(f2filednum,'yyyymmddHHMMSS'),ftell(fid)-(double(frameheaderlen)*4)-(double(extrasamples)*4));
                        index = [index; f2filednum ftell(fid)-(double(frameheaderlen)*4)-(double(extrasamples)*4)];

                        % fast forward in the file by xxx million words (32 bit)
                        % needs to be adjusted depending on sector size!!
                        fseek(fid,200000000*4,0);
                        % find the next frameheader and read it
                        frameheader(2).framesync = fread(fid,1,'*uint32',0,'b');
                        while frameheader(2).framesync ~= FRAMESYNC
                            frameheader(2).framesync = fread(fid,1,'*uint32',0,'b');
                        end

                        frameheader(2).frtimesec = fread(fid,1,'*uint32',0,'b');
                        frameheader(2).frtimeusec = fread(fid,1,'*uint32',0,'b');
                        frameheader(2).lastsmpls = fread(fid,1,'*uint32',0,'b');
                        frameheader(2).turning = fread(fid,1,'*uint32',0,'b');
                        frameheader(2).reserved1 = fread(fid,1,'*uint32',0,'b');
                        frameheader(2).reserved2 = fread(fid,1,'*uint32',0,'b');
                        frameheader(2).reserved3 = fread(fid,1,'*uint32',0,'b');
                        extrasamples = 0;
                    end
                end
            otherwise
                disp('bad state')
        end
        % save frame 2 header as frame 1 for next pass
        frameheader(1) = frameheader(2);
    catch
        % loop craps out at EOF
        fclose(fid);
%         fclose(index_fid);
        save([datapath filename(1:6) '_index.mat'],'index');
        save([datapath filename(1:6) '_info.mat'],'fileheader','minf1turning','maxf1turning','minf2turning','maxf2turning');
        return
    end
end
   