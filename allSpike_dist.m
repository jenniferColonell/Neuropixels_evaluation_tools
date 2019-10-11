function allSpike_dist( varargin )

% input params: fileDir of peaks.mat, fileName of peaks.mat, zMin, zMax, exChan (array)
 
% Z range for histogram, in um:
% set to +/- inf to include all channels
% NP1.0 has vertical pitch = 20 um, 1 bank = 192*20 = 3840 um
% NP2.0 has vertical pitch = 15 um, 1 bank = 192*15 = 2880 um

% exChan:
% Channels to exclude in max, min, average calculations. Should always exclude
% the reference channel (191 for NP 1.0, 127 for NP 2.0.
% Also exclude any known "dead" or "extra noisy" channels on the probe.
% these are given in orignal (zero-based) format

if (length(varargin) == 0)
    % if running standalone, set params here
    zMin = -inf;
    zMax = inf;
    exChan = [127];
    % get thresholding results file from user
    [fileName,fileDir]=uigetfile('*.mat', 'Select threshold results file' );
else
    % params passed in through varargin
    inputCell = varargin(1);
    fileDir = inputCell{1};
    inputCell = varargin(2);
    fileName = inputCell{1};
    inputCell = varargin(3);
    zMin = inputCell{1};
    inputCell = varargin(4);
    zMax = inputCell{1};
    inputCell = varargin(5);
    exChan = inputCell{1};
end

% Reads results *.mat file from detect_merge_peaks.
% Also requires a simple text file of the site coordinates.

% Rebins histograms of spike amplitudes to a common bin size,
% hardcoded to be the binsize of NP data with gain = 500.
% From those histograms, generates mode, 90th, 99th percentile
% for the peak amplitude distribution for each channel.
% That data is output, along with mean threshold, XZ and percCorrect
% for each channel in the file <input name>_chan.txt

% Output a text file of the histogram over select channels: <input name>_hist.txt

    
exChan = exChan + 1;  %conver to 1 based for MATLAB

largeThresh = 200; %threhold for "large events", in uV


% Histogram the raw peaks (res.ampHist) or the merged (res.mAmpHist)
bMerge = 1;


suffPos = strfind(fileName,'_peaks.mat');
outName = [fileName(1:suffPos-1), '_sum.txt'];
outID = fopen( fullfile(fileDir,outName), 'w');
chanOutName = [fileName(1:suffPos-1), '_chan.txt'];
chanOutID = fopen( fullfile(fileDir,chanOutName), 'w');
histName = [fileName(1:suffPos-1), '_hist.txt'];
histID = fopen( fullfile(fileDir,histName), 'w');

load(fullfile(fileDir,fileName));
dataType = res.dataType;

% amplitude histograms are stored in nchanxnbin array res.ampHist
[nChan, nBin] = size(res.ampHist);


% Get site coords from thresholding results file
xPos = res.xPos;
zPos = res.zPos;
shank = res.shank;

inZRange = (zPos > zMin) & (zPos < zMax);

if bMerge
    origHist = res.mAmpHist;
else
    origHist = res.ampHist;
end


removeChan = [];
uVPerBin = 2.34375;      %only exact for neuropixels; recalc if rebinning
switch dataType
    
    case 0
        nBin = 1024;    %amplitudes will be histogrammed over all bits
        uVPerBit = 2.34375;  
        sumBin = 1;
        uVPerBin = sumBin*uVPerBit;
        goodChan = 0;    
        for i = 1:nChan
            isBad = sum(find(exChan == i));            
            if( inZRange(i) == 1 && isBad == 0 )
                goodChan = goodChan + 1;
                tempHist(goodChan,:) = origHist(i,:);
                newZ(goodChan) = zPos(i);
                newX(goodChan) = xPos(i);
                newShank(goodChan) = shank(i);
            else
                removeChan = [removeChan,i];
            end
        end         
        res.ampHist = [];
        res.ampHist = tempHist;
        size(res.ampHist);
        xPos = [];
        xPos = newX;
        zPos = [];
        zPos = newZ;
        shank = [];
        shank = newShank;

    case 1
        nBin = 16384;    %amplitudes will be histogrammed over all bits
        uVPerBit = 0.7629;
        goodChan = 0;
        for i = 1:nChan
            isBad = sum(find(exChan == i));
            if( inZRange(i) == 1 && isBad == 0 )
                goodChan = goodChan + 1;
                tempHist(goodChan,:) = origHist(i,:);
                newZ(goodChan) = zPos(i);
                newX(goodChan) = xPos(i);
                newShank(goodChan) = shank(i);
            else
                removeChan = [removeChan,i]; 
            end
        end  
        
        res.ampHist = [];
        res.ampHist = tempHist;
        %size(res.ampHist)
        xPos = [];
        xPos = newX;
        zPos = [];
        zPos = newZ;
        shank = [];
        shank = newShank;
        
        %now rebin the histograms so the bins are as close as possible to
        %NP1.0 bins        
        sumBin = round(uVPerBin/uVPerBit);
        uVPerBin = sumBin*uVPerBit;
        nNewBin = floor(nBin/sumBin);
        tempHist = zeros(goodChan,nNewBin,'int32');
        %size(tempHist)
        for k = 1:nNewBin
            currStart = (k-1)*sumBin + 1;
            for i = 1:goodChan
                tempHist(i,k) = sum(res.ampHist(i,(currStart:currStart+sumBin-1)));
            end
        end      
        nBin = nNewBin;
        res.ampHist = [];
        res.ampHist = tempHist;
        
    otherwise
        fprintf( 'unknown dataType\n' );
        return;
end

MAD_goodChan = res.MAD;
MAD_goodChan(removeChan) = [];
goodChanList = (0:383);
goodChanList(removeChan) = [];

pVec = [0.5, 0.9, 0.99];
xVal = uVPerBin*(0:nBin-1);

%loop over channels, write out threshold, nSpike, mode of amplitudes, 90th perc
%accumulate number of channels with 99th > "largeThresh"
nLarge = 0;

fprintf( chanOutID, 'shank\tX\tZ\t' );
fprintf( chanOutID, 'chan\test RMS\teventRate\tampMod\t90th ptile\t99th ptile\n');

currRate = zeros(size(MAD_goodChan));
curr90th = zeros(size(MAD_goodChan));
curr99th = zeros(size(MAD_goodChan));

for i = 1:goodChan
    [maxCnt, maxI] = max(res.ampHist(i,:));
    currMode = xVal(maxI);
    curr90th(i) = xVal(pctileFromHist(res.ampHist(i,:),0.9));
    curr99th(i) = xVal(pctileFromHist(res.ampHist(i,:),0.99));    
    currRate(i) = sum(res.ampHist(i,:))/res.analyzedSec;
    fprintf( chanOutID, '%d\t%d\t%d\t', shank(i), xPos(i), zPos(i) );
    fprintf(chanOutID,'%d\t%.3f\t%.3f\t%.3f\t',goodChanList(i),MAD_goodChan(i),currRate(i),currMode);
    fprintf(chanOutID,'%.3f\t%.3f\n',curr90th(i),curr99th(i));
    if curr99th(i) > largeThresh
        nLarge = nLarge + 1;
    end
end
fclose(chanOutID);

%for plotting resort in z, then x
[~,zsorder] = sortrows([zPos;shank]',[2,1]);

% checking plot order
%  figure(101)
%  plot(1:goodChan, zPos(zsorder))
%  figure(102)
%  plot(1:goodChan, xPos(zsorder))

%sum along channel to get full histogram
allChan = sum(res.ampHist,1)';    
totalCounts = sum(allChan);

normHist = allChan/totalCounts;



f1 = figure('Units','Centimeters','Position',[3,12,8,8]);
plot(xVal, allChan);
titleStr = sprintf('Amplitude Distribution for %s', fileName );
title(titleStr, 'Interpreter', 'none');
xlabel('amplitude (uV)');
ylabel('number of spikes');

% bulid inverse CDF map for each channel
% code from Nick Steinmetz,
f1 = figure('Units','Centimeters','Position',[4,11,8,8]);
cdf = single(cumsum(fliplr(res.mAmpHist),2));
cdf = cdf/res.analyzedSec;  %convert to spike rates
[nRow, nCol] = size(cdf);

%flip back
cdf = fliplr(cdf);

ampLimit = 1000;
for i = 1:numel(xVal)
    if xVal(i) > ampLimit
        xMax = i;
        break;
    end
end

depthX = 1:nRow;
imagesc(xVal(1:xMax), depthX, cdf(:,(1:xMax)));
xlabel('spike amplitude (uV)');
ylabel('channel index');
title('inverse cdf');
set(gca,'YDir','normal');
colorbar
colormap(colormap_greyZero_blackred)
caxis([0 20]);

f3 = figure('Units','Centimeters','Position',[8,5,14,18]);
subplot(3,1,1);
plot(1:numel(goodChanList), MAD_goodChan(zsorder));
titleStr = sprintf('est RMS' );
title(titleStr, 'Interpreter', 'none');
xlabel('channels sorted by shank, then z');
ylabel('estimated rms (uV)');

subplot(3,1,2);
plot(1:numel(goodChanList), currRate(zsorder));
titleStr = sprintf('Event Rate' );
title(titleStr, 'Interpreter', 'none');
xlabel('channels sorted by shank, then z');
ylabel('event rate (Hz)');

subplot(3,1,3);
plot(1:numel(goodChanList), curr90th(zsorder), 1:numel(goodChanList), curr99th(zsorder));
titleStr = sprintf('90th (blue) and 99th (orange) percentile' );
title(titleStr, 'Interpreter', 'none');
xlabel('channels sorted by shank, then z');
ylabel('voltage at percentile (uV)');


if isfield(res, 'probeType' )
    plotChanFiles( res.probeType, fullfile(fileDir,chanOutName) );
end

fprintf( histID, 'uV\thistCounts\tnormHist\n');
for i = 1:nBin
    fprintf(histID, '%.3f\t%d\t%.6e\n', xVal(i), allChan(i), normHist(i));
end

% Calculate summary values
eventRate = sum(allChan)/res.analyzedSec;

avgMAD = mean(MAD_goodChan);
stdMAD = std(MAD_goodChan);

% get percentile values
for i = 1:numel(pVec)
    pctInd(i) = pctileFromHist( allChan, pVec(i) );
end

% channels with 99th percentile events > 200 uV (how many big units do we
% see?)


fprintf( outID, 'fileName\teventRate\t50th\t90th\t99th\t#channels with 99th > 200\test rms\tstd\tmax\tmin\n' );
fprintf( outID, '%s\t%.3f\t%.3f\t%.3f\t%.3f\t%d\t%.3f\t%.3f\t%.3f\t%.3f\n', ...
      fileName(1:suffPos-1), eventRate, xVal(pctInd(1)), xVal(pctInd(2)), xVal(pctInd(3)), ...
      nLarge, avgMAD, stdMAD, max(MAD_goodChan), min(MAD_goodChan) );

%fprintf(outID, '\n\n\n' );

%Next, get number of channels that contribute to spikes above a given value
%Need to translate the pctInd values back to original indicies in the full
%(un-rebinned) histograms

% fprintf( outID, 'percent\tnContrib\n');
% for i = 1:numel(pVec)
%     [contribChans, cProfile] = chanProfile( res.ampHist, pctInd(i) );
%     nContrib = numel(contribChans);
%     fprintf( outID, '%.1f\t%d', 100*pVec(i), nContrib );
%     if( numel(contribChans) < 500 )
%         [cProfile,sortI] = sort(cProfile,'descend');
%         contribChans = contribChans(sortI);
%         for j = 1:nContrib
%             fprintf( outID, '\t%d', contribChans(j) );
%         end
%         fprintf( outID, '\n\t');
%         for j = 1:nContrib
%             fprintf( outID, '\t%.6e', cProfile(j) );
%         end
%     end
%     fprintf( outID, '\n');
% end

fclose(histID);
fclose(outID);

end

function [histPctile] = pctileFromHist( hData, p )

%calculate the bin corresponding to a given percentile p (value 0 to 1)
%from a histogram

nBin = numel(hData);

if sum(hData) > 0
    targVal = p*sum(hData);

    i = 1;
    currVal = 0;

    while( currVal <= targVal && i < nBin )
        currVal = currVal + hData(i);
        i = i + 1;
    end

    histPctile = int32(i);
else
    histPctile = 1;
end
  
end

function [chans, cProfile] = chanProfile( ampHist, pctInd )
    %ampHist is [nChan x nBin]
    %pctInd is a bin index
    %get chan indicies and contributions for all spikes >= bin Index
    [nChan, nBin] = size(ampHist);
    totHist = sum(ampHist,1)';
    pctInd = int32(pctInd);
    totalArea = sum(totHist(pctInd-1:nBin));
    nContrib = 0;
    chans = [];
    cProfile = [];
    for i = 1:nChan
        currSum = sum(ampHist(i,pctInd-1:nBin));
        if currSum > 0
            nContrib = nContrib + 1;
            chans(nContrib) = i;
            cProfile(nContrib) = currSum/totalArea;
        end
    end
end

function map = colormap_greyZero_blackred()
    % custom colormap for inverse cdf plot from Nick Steinmetz
    map = hot(1000);
    map(1,:) = [0.4,0.4,0.4];
end