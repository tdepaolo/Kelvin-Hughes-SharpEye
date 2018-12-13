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
% [sdnum, slat, slon, scog, ssog, shead] = load_shipData(datapath(4:11));
% MakeQTMovie quality 1
% MakeQTMovie framerate 1
% % MakeQTMovie size [6 6]
% MakeQTMovie start SharpEyeScans.qt

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
     lonp,latp] = read_SharpEye_sector([datapath datafile], index(ns,2));
    toc
    
    if ploton
        figure(1)
%         subplot(1,2,1)
        pcolor(lonp,latp,(nanmean(double(f1mag),3)));
        shading flat
        axis square
        colorbar
        caxis([0 8])
        title(datestr(scandnum,'yyyymmddHHMMSS'));
        drawnow;
        
        figure(2)
%         subplot(1,2,1)
        pcolor(lonp,latp,(nanmean(double(f2mag),3)));
        shading flat
        axis square
        colorbar
        caxis([0 8])
        title(datestr(scandnum,'yyyymmddHHMMSS'));
        drawnow;
%         saveas(gcf,[datapath 'scan_' datestr(scandnum,'yyyymmddHHMMSS') '.jpg'],'jpg');
%         MakeQTMovie addfigure
%         close(1);
    end
end
% MakeQTMovie finish
    
    
      