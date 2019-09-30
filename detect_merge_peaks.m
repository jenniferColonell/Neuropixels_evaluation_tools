function detect_merge_peaks()

% Set the dataType and qqFactor before running
% 0 for neuropixels 1.0, 1 for NP 2.0
dataType = 0;

%thresholding params
bUseConstThresh = 1;
constThresh = 80;    %in uV
qqFactor = 5; % for JRClust-style thresholding; ignored if bUseConstThresh = 1

%other params for JRClust type merging
hCfg.evtDetectRad = 50; %in um, distance to look for duplicate peaks
hCfg.refracIntSamp = 7; %in samples, distance in time to look for duplicat peaks

% for either data type
nchan = 385;
dataChan = 384;

maxTime = 600; %in sec; takes the first maxTime secons in the file

% get data file from user. 
[fileName,fileDir]=uigetfile('*.bin', 'Select file' );

cd(fileDir);

% Build name for output file, which will be a matlab structure
[~,inName,~] = fileparts(fileName);
outName = [inName,'_peaks.mat'];



% get coordinate file from user, chan, X, Y, tab delimited
% X coordiantes across shanks should reflect the distance between the
% shanks. should include all channels in teh file (e.g. include the ref
% channels, if present).
% Get site coords file from user
[coordName,coordDir]=uigetfile('*.txt', 'Select coords file' );
cID = fopen( fullfile(coordDir,coordName), 'r' );
hCfg.siteLoc = zeros(dataChan,2);
for i = 1:dataChan
    tline = fgetl(cID);
    currDat = sscanf(tline, '%d%d%d');
    hCfg.siteLoc(i,1) = currDat(2);
    hCfg.siteLoc(i,2) = currDat(3);
end
fclose(cID);




switch dataType
    
    case 0
        fs = 30000;
        nBit = 1024;    %amplitudes will be histogrammed over all bits
        uVPerBit = 2.34375;
     
    case 1
        fs = 30000;
        nBit = 16384;    %amplitudes will be histogrammed over all bits
        uVPerBit = 0.7629;
        
    otherwise
        fprintf( 'unknown dataType\n' );
        return;
end

cThreshBits = constThresh/uVPerBit;

edges = (0:nBit);   %nBit+1 edges -> nBit bins. Overflow goes to top/bottom bin


    fileStats = dir(fullfile(fileDir,fileName));

    fileSize = fileStats.bytes;
    
    fm = memmapfile(fullfile(fileDir,fileName),'Format','int16');
    
    maxBatchBytes = 1e8;
    maxBatchSamples = floor(maxBatchBytes/(nchan*2));
    batchBytes = maxBatchSamples*nchan*2;
    
    nneighBelow = 1;

    runSec = fileSize/(nchan*2*fs);
    fprintf('Run length in seconds: %.2f\n',runSec);

    fileSamples = fileSize/(nchan*2);
    
    if runSec < maxTime
        nBatch = floor(fileSamples/maxBatchSamples) + 1;
    else
        useSamples = maxTime*fs;
        nBatch = floor(useSamples/maxBatchSamples) + 1;
    end
        
    batchSamples = maxBatchSamples;

    fprintf('nBatch, batchSamples: %d, %d\n', nBatch, batchSamples );
    
    %will be processing the nBatch-1 "full" batches (saves complexity)
    analyzedSec = (nBatch-1)*batchSamples/fs;

    sizeA = [nchan,batchSamples];

    %set up arrays to hold amplitude histograms for each channel
    %for NP1.0 data, max amplitude is 10 bits; for NP 2.0, max is 14 bits

    ampHist = zeros(dataChan,nBit,'int32');   
    tempHist = zeros(1,nBit,'int32');
    
    %Note -- should we try to preallocate space for these?
    %Would need a way to estimate the number of spikes
    %Could consider stopping after some number of spikes are reached.
    allTimes = [];
    allAmps = [];
    allSites = [];

    %likewise save thresholds
    siteThresh = zeros(dataChan, nBatch-1, 'single'); %we'll be ignoring the final partial batch

    % 
    % For each batch:
    %     read data
    %     loop over sites
    %         calculate threshold
    %         find peaks
    %         get bipolar peak amplitude
    %         add to histogram of amplitudes for that site
    %         
    %     get median values for each channel
    % 
    %   

    % for simplicity, just reading in whole batches, not using the last 
    % segment of data.
    % also, assuming the data is at least as large as one batch.
    
tic
    batch16bit = batchBytes/2;
    for i = 1:nBatch-1
         fprintf( 'Reading batch: %d of %d\n', i, nBatch-1 );
         batchOffset = int64((i-1)*batch16bit) + 1;
         dataArray = fm.Data(batchOffset:batchOffset+batch16bit-1);
         dataArray = reshape(dataArray,sizeA);
        
         % Calculate thresholds for this batch, want to get the MAD result
         % even if using constant threshold
         for jChan = 1:dataChan
            siteThresh(jChan,i) = calcThresh(dataArray(jChan,:), qqFactor, nBit);
         end

         %for each channel, find peaks (using JRClust algorithm) and add into
         %histogram   
         
         batchPeaks = 0;
         batchStart = (i-1)*batchSamples + 1; %one based for matlab
         for jChan = 1:dataChan
             if bUseConstThresh
                thresh = cThreshBits;
             else
                thresh = siteThresh(jChan,i);
             end
             [peaks, amps] = findPeaks(dataArray(jChan,:), thresh, nneighBelow, fs);
             allAmps = [allAmps,amps];
             allTimes = [allTimes,peaks+batchStart];            
             currSites = jChan*ones(1,numel(peaks));
             allSites = [allSites,currSites];
             [tempHist, ~] = (histcounts(amps, edges));
             ampHist(jChan,:) = ampHist(jChan,:) + int32(tempHist);
             batchPeaks = batchPeaks + numel(peaks);

         end
         fprintf('number of peaks: %d\n', batchPeaks);
         clear('dataArray');
    end

    
fprintf( "Time to run detection: %.3f\n", toc);
fprintf( "total number of spikes: %d\n", numel(allTimes));

    %save result structure
    res.cThresh = constThresh;
    res.useConstThresh = bUseConstThresh;
    res.Thresh = siteThresh;
    res.ampHist = ampHist;
    res.qqFactor = qqFactor;
    res.analyzedSec = analyzedSec;
    
    %count up number of peaks and calculate mean MAD for each channel
    for i = 1:dataChan
        res.nPeak(i) = sum(ampHist(i,:));
        res.MAD(i) = uVPerBit*mean(res.Thresh(i,:))/qqFactor;
    end
    
    %start a parallel pool to run the merge step
    %delete any currently running pool
    %delete(gcp('nocreate'))
    %create parallel pool
    if isempty(gcp)
    locCluster = parcluster('local');
    parpool('local',locCluster.NumWorkers);
    end
    [spikeTimes, spikeAmps, spikeSites] = mergePeaks(allTimes, allAmps, allSites, dataChan, hCfg);
    
    res.mergedTimes = spikeTimes;
    res.mergedAmps = spikeAmps;
    res.mergedSites = spikeSites;
    %build histograms of merged peaks
    mAmpHist = zeros(dataChan,nBit,'int32');
    for jChan = 1:dataChan
        amps = spikeAmps(spikeSites == jChan);
        [tempHist, ~] = (histcounts(amps, edges));
        mAmpHist(jChan,:) = mAmpHist(jChan,:) + int32(tempHist);
    end
    res.mAmpHist = mAmpHist;
    save( fullfile(fileDir,outName), 'res' );

end


function thresh = calcThresh( currSamples, qqFactor, nBit )

    quirogaDenom = 0.6745;
    
    currDev = abs(currSamples);
    
    nEdge = nBit;
    hist_edges = 0:nEdge;
    [counts,~] = histcounts(currDev,hist_edges);
    nValues = sum(counts);
    
    medFound = 0;
    iE = 1;
    while( ~medFound && iE < nEdge )
        sumLow = sum(counts(1:iE));
        sumHigh = sum(counts(iE+1:nEdge));
        if( sumHigh < sumLow )
            medFound = 1;
            currMed = hist_edges(iE);
            %calculate the estimated median within the median bin using a
            %linear interpolation between the lower and upper edge of the
            %median bin. This is known as estimating the mean for grouped
            %data -- the groups in this case are the low bits of the int16
            %values.
            
            lb = currMed; %lower bound of bin containing median
            if iE > 1
                yDist = nValues/2 - sum(counts(1:iE-1));
            else
                yDist = nValues/2; %if median is in zero bin, sum of bins below med = 0;
            end
            estMed = lb + yDist/counts(iE);
            
            %fprintf( '%d\t%d\t%.2f\n', selectChan(i), currMed, estMed(selectChan(i)) );
        else
            iE = iE + 1;
        end
    end
    
    if( medFound == 0 ) 
            fprintf( 'median not found. problem reading data?\n' );
            thresh = -1;
            return;
    else
        thresh = qqFactor*(estMed/quirogaDenom);
    end

end



% from JRclust 4.0.0-alpha.5
function [peaks, amps] = findPeaks(samplesIn, thresh, nneighBelow, fs)
    %FINDPEAKS Find samples which exceed threshold
    %params for calculating amplitudes
    sT = int32(floor((0.25/1000)*fs));  % req time after batch start
    eT = int32(floor((0.4/1000)*fs));    % req time before batch end
    %amplitude is max - min over the span -st to +eT
    
    if isempty(nneighBelow)
        nneighBelow = 1;
    end

    peaks = [];
    amps = [];
    if isempty(samplesIn)
        return;
    end

    %exceedsThresh = jrclust.utils.tryGather(samplesIn < -abs(thresh));
    exceedsThresh = (samplesIn < -abs(thresh));
    peakLocs = find(exceedsThresh);
    if isempty(peakLocs)
        return;
    end

    % allow only peaks far enough from start to get complete waveform
    ind = find(peakLocs < sT + 1);
    if numel(ind) == numel(peakLocs) % these are the only peaks found, none valid
        return;
    else % discard all with times less than sT + 1
        peakLocs(ind) = [];
    end

    % allow only peaks far enough from end to get complete waveform
    ind = find(peakLocs > numel(samplesIn) - eT );
    if numel(ind) == numel(peakLocs) % these are the only peaks found, none valid
        return;
    else % discard all with times greater than eT
        peakLocs(ind) = [];
    end

    peakCenters = samplesIn(peakLocs);
    % take only "peaks" whose sample neighbors are not larger
    peaks = peakLocs(peakCenters <= samplesIn(peakLocs + 1) & peakCenters <= samplesIn(peakLocs - 1));
    if isempty(peaks)
        return;
    end

    % take only peaks who have one or two sample neighbors exceeding threshold
    if nneighBelow == 1
        peaks = peaks(exceedsThresh(peaks - 1) | exceedsThresh(peaks + 1));
    elseif nneighBelow == 2
        peaks = peaks(exceedsThresh(peaks - 1) & exceedsThresh(peaks + 1));
    end
    
    %single point amplitudes
    amps = zeros(1,numel(peaks));
    for i = 1:numel(peaks)
        amps(i) = abs(samplesIn(peaks(i)));
    end
    
% alternate version: max - min over nearest 0.5 msec
%     for i = 1:numel(peaks)
%         amps(i) = max(samplesIn(peaks(i)-sT:peaks(i)+eT)) - min(samplesIn(peaks(i)-sT:peaks(i)+eT));
%     end
end


function [spikeTimes, spikeAmps, spikeSites] = mergePeaks(allTimes, allAmps, allSites, nSites, hCfg)
    %adapted from JRC, alpha 4.0, commit 4280a02
    %MERGEPEAKS Merge duplicate peak events
    spikeTimes = allTimes';
    spikeAmps = allAmps';
    spikeSites = allSites';

    [spikeTimes, argsort] = sort(spikeTimes);
    spikeAmps = spikeAmps(argsort);
    spikeSites = int32(spikeSites(argsort));
    spikeTimes = int32(spikeTimes);

    [mergedTimes, mergedAmps, mergedSites] = deal(cell(nSites,1));

    try
        parfor iSite = 1:nSites
            try
                [mergedTimes{iSite}, mergedAmps{iSite}, mergedSites{iSite}] = ...
                    mergeSpikesSite(spikeTimes, spikeAmps, spikeSites, iSite, hCfg);
            catch ME% don't try to display an error here
               warning('failed to use parallel pool %d: %s', iSite, ME.message);
            end
        end
    catch % parfor failure
        for iSite = 1:nSites
            try
                [mergedTimes{iSite}, mergedAmps{iSite}, mergedSites{iSite}] = ...
                    mergeSpikesSite(spikeTimes, spikeAmps, spikeSites, iSite, hCfg);
            catch ME
                warning('failed to merge spikes on site %d: %s', iSite, ME.message);
            end
        end
    end

    % merge parfor output and sort
    spikeTimes = jrclust.utils.neCell2mat(mergedTimes);
    spikeAmps = jrclust.utils.neCell2mat(mergedAmps);
    spikeSites = jrclust.utils.neCell2mat(mergedSites);

    [spikeTimes, argsort] = sort(spikeTimes); % sort by time
    spikeAmps = jrclust.utils.tryGather(spikeAmps(argsort));
    spikeSites = spikeSites(argsort);
end

%% LOCAL FUNCTIONS
function [timesOut, ampsOut, sitesOut] = mergeSpikesSite(spikeTimes, spikeAmps, spikeSites, iSite, hCfg)
    %MERGESPIKESSITE Merge spikes in the refractory period
    nLims = int32(abs(hCfg.refracIntSamp));

    % find neighboring spikes
    nearbySites = jrclust.utils.findNearbySites(hCfg.siteLoc, iSite, hCfg.evtDetectRad); % includes iSite
    spikesBySite = arrayfun(@(jSite) find(spikeSites == jSite), nearbySites, 'UniformOutput', 0);
    timesBySite = arrayfun(@(jSite) spikeTimes(spikesBySite{jSite}), 1:numel(nearbySites), 'UniformOutput', 0);
    ampsBySite = arrayfun(@(jSite) spikeAmps(spikesBySite{jSite}), 1:numel(nearbySites), 'UniformOutput', 0);

    iiSite = (nearbySites == iSite);
    iSpikes = spikesBySite{iiSite};
    iTimes = timesBySite{iiSite};
    iAmps = ampsBySite{iiSite};

    % search over peaks on neighboring sites and in refractory period to
    % see which peaks on this site to keep
    keepMe = true(size(iSpikes));
    for jjSite = 1:numel(spikesBySite)
        jSite = nearbySites(jjSite);

        jSpikes = spikesBySite{jjSite};
        jTimes = timesBySite{jjSite};
        jAmps = ampsBySite{jjSite};

        if iSite == jSite
            delays = [-nLims:-1, 1:nLims]; % skip 0 delay
        else
            delays = -nLims:nLims;
        end

        for iiDelay = 1:numel(delays)
            iDelay = delays(iiDelay);

            [jLocs, iLocs] = ismember(jTimes, iTimes + iDelay);
            if ~any(jLocs)
                continue;
            end

            % jLocs: index into j{Spikes,Times,Amps}
            % iLocs: (1st) corresponding index into i{Spikes,Times,Amps}/keepMe
            jLocs = find(jLocs);
            iLocs = iLocs(jLocs);

            % drop spikes on iSite where spikes on jSite have larger
            % magnitudes
            nearbyLarger = abs(jAmps(jLocs)) > abs(iAmps(iLocs));
            keepMe(iLocs(nearbyLarger)) = 0;

            % flag equal-amplitude nearby spikes
            ampsEqual = (jAmps(jLocs) == iAmps(iLocs));
            if any(ampsEqual)
                if iDelay < 0 % drop spikes on iSite occurring later than those on jSite
                    keepMe(iLocs(ampsEqual)) = 0;
                elseif iDelay == 0 && jSite < iSite % drop spikes on iSite if jSite is of lower index
                    keepMe(iLocs(ampsEqual)) = 0;
                end
            end
        end
    end

    % keep the peak spikes only
    timesOut = iTimes(keepMe);
    ampsOut = iAmps(keepMe);
    sitesOut = repmat(int32(iSite), size(timesOut));
end
