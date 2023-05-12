function meas_peak_to_peak( varargin )
% Run on UNFILTERED data, to get median, max, min, and peak to peak 
% voltage. peak_to_peak measurements characterize the dynamic range used; 
% the median is the voltage offset. Primary output is file 
% <run_name>_PtoP_chan.txt (also a mat file of the same data).
% Call SGLXMetaToCoords first to make coordinate file.

if (length(varargin) == 0)
    % probe type: NP1.0 = 1; NP2.0 SS = 21, NP2.0 = 24;
    probeType = 2013;
else
    inputCell = varargin(1);
    probeType = inputCell{1};
end

% set dataType
% 0 for neuropixels 1.0, 1 for NP 2.0. 
if probeType == 1
    dataType = 0;
elseif (probeType == 21) || (probeType == 24)
    dataType = 1;
elseif probeType == 2013
    dataType = 2;
else
    fprintf( "unknown probe type\n" );
    return
end

% number of data channels; if saving all, nchan = 385, dataChan = 384
nchan = 385;
dataChan = 384;

maxTime = 10; %in sec; takes the first maxTime seconds in the file

%meas_params
avgTime = 1; %in sec; this will be the batch size
largeV = 150; %in uV, channels with p-p larger than this are counted as large

% bUseConstThresh = 1;
% constThresh = 80;    %in uV
% qqFactor = 5; % for JRClust-style thresholding; ignored if bUseConstThresh = 1

% %other params for JRClust type merging
% hCfg.evtDetectRad = 50; %in um, distance to look for duplicate peaks
% hCfg.refracIntSamp = 7; %in samples, distance in time to look for duplicat peaks

% get data file from user. 
[fileName,fileDir]=uigetfile('*.bin', 'Select file' );

cd(fileDir);

% Build name for output file, which will be a matlab structure
[~,inName,~] = fileparts(fileName);
outName = [inName,'_PtoP.mat'];



% get coordinate file from user, chan, X, Y, shank index, tab delimited
% X coordiantes across shanks should reflect the distance between the
% shanks. should include all channels in teh file (e.g. include the ref
% channels, if present).
% Get site coords file from user
[coordName,coordDir]=uigetfile('*.txt', 'Select coords file' );
cID = fopen( fullfile(coordDir,coordName), 'r' );
hCfg.siteLoc = zeros(dataChan,2);
shank = zeros(dataChan,1);
for i = 1:dataChan
    tline = fgetl(cID);
    currDat = sscanf(tline, '%d%d%d%d');
    hCfg.siteLoc(i,1) = currDat(2);
    hCfg.siteLoc(i,2) = currDat(3);
    shank(i) = currDat(4);
end
fclose(cID);


switch dataType
    
    case 0
        fs = 30000;
        nBit = 1024;    %amplitudes will be histogrammed over all bits
        uVPerBit = 2.34375; %assumes gain is set to 500
     
    case 1
        fs = 30000;
        nBit = 16384;    %amplitudes will be histogrammed over all bits
        uVPerBit = 0.7629;

    case 2
        fs = 30000;
        nBit = 4096;
        uVPerBit = 1e6*((2*0.62)/100)/(2*2048); % range = +/- 0.62 V, gain = 100, maxInt = 2048 => 2*2048 bits total

    otherwise
        fprintf( 'unknown dataType\n' );
        return;
end




    fileStats = dir(fullfile(fileDir,fileName));

    fileSize = fileStats.bytes;
    
    fm = memmapfile(fullfile(fileDir,fileName),'Format','int16');
    
    % set batch size to avgTime set by user. Make sure the setting is 
    % not absurdly large.
    maxBatchBytes = 1e8;
    maxBatchSamples = floor(maxBatchBytes/(nchan*2));
    
    batchSamples = avgTime*fs;
    batchBytes = batchSamples*nchan*2;
    
    if batchBytes > maxBatchBytes
        fprintf( 'avgTime of %.1f is too large. Set a smaller avgTime, e.g. 1 second\n', avgTime );
        return;
    end
      
    
    runSec = fileSize/(nchan*2*fs);
    fprintf('Run length in seconds: %.2f\n',runSec);

    fileSamples = fileSize/(nchan*2);
    
    if runSec < maxTime
        nBatch = floor(fileSamples/batchSamples) + 1;
    else
        useSamples = maxTime*fs;
        nBatch = floor(useSamples/batchSamples) + 1;
    end
        
    fprintf('nBatch, batchSamples: %d, %d\n', nBatch, batchSamples );
    
    %will be processing the nBatch-1 "full" batches (saves complexity)
    analyzedSec = (nBatch-1)*batchSamples/fs;

    sizeA = [nchan,batchSamples];

    % make arrays to hold max, min, peak-peak, median for each batch

    maxV = zeros(dataChan, nBatch-1, 'single');
    minV = zeros(dataChan, nBatch-1, 'single');
    ppV = zeros(dataChan, nBatch-1, 'single');
    medV = zeros(dataChan, nBatch-1, 'single'); 
    
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
            %medV(jChan,i) = uVPerBit*calcThresh(dataArray(jChan,:), 1, nBit);
            medV(jChan,i) = uVPerBit*median(dataArray(jChan,:));
         end
         maxV(:,i) = uVPerBit*single(max(dataArray(1:dataChan,:),[],2));
         minV(:,i) = uVPerBit*single(min(dataArray(1:dataChan,:),[],2));
         ppV(:,i) = maxV(:,i) - minV(:,i);
         clear('dataArray');
    end

    
fprintf( "Time to run measurement: %.3f\n", toc);

    %save result structure
    res.maxV = maxV;
    res.minV = minV;
    res.medV = medV;
    res.ppV = ppV;
    res.analyzedSec = analyzedSec;
    res.dataType = dataType;
    res.probeType = probeType;
    res.xPos = hCfg.siteLoc(:,1);
    res.zPos = hCfg.siteLoc(:,2);
    res.shank = shank;

   
    % get averages. Could do this over the array, put in loop in case
    % later want to exclude outliers
    for i = 1:dataChan
        res.avg_ppV(i) = mean(ppV(i,:));
        res.std_ppV(i) = std(ppV(i,:));
        res.MAD(i) = mean(medV(i,:));
    end
    save( fullfile(fileDir,outName), 'res' );
    
    % write out a chan file that will be readable by plotChanFiles
    chanOutName = [inName, '_PtoP_chan.txt'];
    chanOutID = fopen( fullfile(fileDir,chanOutName), 'w');
    fprintf( chanOutID, 'shank\tX\tZ\t' );
    fprintf( chanOutID, 'chan\tmedian\tPtoP\tPtoP std\n');
    for i = 1:dataChan
        fprintf( chanOutID, '%d\t%d\t%d\t', shank(i), hCfg.siteLoc(i,1), hCfg.siteLoc(i,2) );
        fprintf( chanOutID,'%d\t%.3f\t%.3f\t%.3f\n',i,res.MAD(i),res.avg_ppV(i),res.std_ppV(i));
    end
    fclose(chanOutID);
    
    %call probe plotter (commented out because it isn't very useful).
%     plotChanFiles( res.probeType, 1, fullfile(fileDir,chanOutName) );
    
%     sumOutName = [inName, '_PtoP_sum.txt'];
%     sumOutID = fopen( fullfile(fileDir,sumOutName), 'w');
% 
%     nLarge = sum(res.avg_ppV > largeV);
%     fprintf( sumOutID, 'fileName\tmax PtoP\tchan with PtoP > %.0f\n',largeV);
%     fprintf( sumOutID, '%s\t%.1f\t%d\n', inName, max(res.avg_ppV), nLarge); 
%     fclose(sumOutID);

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
            fprintf( 'median not found. problem reading data or all datapoints equal\n' );
            thresh = -1;
            return;
    else
        thresh = qqFactor*(estMed/quirogaDenom);
    end

end



