


%% initialize params
clear
close all

NUM_CHANNELS = 256;
SAMPLERATE = 20000;
totalTimeStim = 20;
timeStim = 10;
reps = 25;

spotSize = 100;

freq = 25;%floor(1000/timeStim); %floor(1000/totalTimeStim);


considerWindowMs = 200; %ms Window to consider post stim
considerWindowInd = considerWindowMs*SAMPLERATE/1000; %num_points Window to consider post stim

time = -considerWindowMs:considerWindowMs; %ms Window to consider post stim

ChannelsToConsider = [3:14 16:127 129:254];


MEA_MAP = [ 127 130 223 254 55 91 122 21 52 88 115 18 45 81 82 128 ;...
            196 226 193 224 25 56 92 119 87 118 17 48 84 111 112 13 ;...
            195 225 158 194 123 26 53 89 117 20 47 83 114 44 79 109 ;...
            159 198 228 157 93 124 23 54 19 50 86 113 43 80 110 11 ;...
            131 162 197 227 253 94 121 24 49 85 116 14 77 107 12 39 ;...
            229 134 161 200 230 160 156 51 22 15 46 78 108 9 40 75 ;...
            202 232 133 164 199 132 155 90 120 41 42 105 10 37 76 106 ;...
            166 201 231 136 163 135 234 129 16 74 35 07 38 73 103 8 ;...
            138 233 203 168 137 165 204 215 125 104 5 33 6 101 71 36 ;...
            236 206 167 140 235 174 173 150 186 63 02 69 34 3 102 72 ;...
            205 170 139 238 208 183 245 185 151 247 30 100 70 31 04 99 ;...
            169 144 237 207 142 180 192 153 59 250 220 217 97 67 32 01 ;...
            143 242 212 171 182 145 175 187 149 57 248 218 27 98 68 29 ;...
            241 211 172 184 243 177 189 152 62 251 221 246 216 28 95 65 ;...
            141 240 239 213 179 191 176 188 148 60 252 222 244 214 96 66 ;...
            255 210 209 181 147 178 190 154 61 146 58 249 219 64 126 256 ]; 

        
% MEA_MAP1 = MEA_MAP1(:, end:-1:1);
% 
% MEA_MAP1 = MEA_MAP;
% % % MEA_MAP1 = MEA_MAP1(end:-1:1, :);
% MEA_MAP1 = MEA_MAP1(:, end:-1:1);
% MEA_MAP1 = MEA_MAP1';
% MEA_MAP1 = MEA_MAP1(:);

% ELEC2MEA_MAP = zeros(1, 256);
% for i=1:256
% 	ELEC2MEA_MAP(i) = electrode(i).ID == MEA_MAP;
% end


sampling_rate = 20000;

%% load vec and frametime file

disp('loading data..');

date_ = '11';

pathBase = '/media/icub/B4A68AB8A68A7B1E/icub/POST_DOC/';     
pathMain = [pathBase '/DATA/' date_ '_10_2016_opto_primate/']; % initialize path of data location

load([pathMain '/' 'cacahuete_day_20161011_EWS__.mat']); % load sorted cell data

load([pathMain '/' 'test0_11-Oct-2016_17h21.mat']); % load workspace containing stimulation parameters

[file, pathDir] = uigetfile([pathMain '/*.raw'],'Select any channel file:' , 'MultiSelect', 'on'); % select the raw file to get spike shapes
indexEnd = strfind(file, '_C_');

%%

f = size(file, 1);

THRESHOLD = 50;

h_r = mcd_header([pathDir file]);

m = memmapfile([pathDir file], 'Offset', h_r.data_start, 'Format', 'uint16');
%map the raw file to the variable m such that m.Data contains the raw
%signal

% h_r = mcd_header([pathDir Filename{f} '.raw']);
% 
% m = memmapfile([pathDir Filename{f} '.raw'], 'Offset', h_r.data_start, 'Format', 'uint16');

Time = h_r.npts/SAMPLERATE;

[b,a] = butter(2, 1000/(SAMPLERATE /2), 'high');

%%
%tests = {'EWS_100um_100hz', 'EWS_50um_100hz', 'EWS_25um_100hz', 'EWS_100um_50hz', 'EWS_50um_50hz', 'EWS_25um_50hz', 'EWS_100um_25hz', 'EWS_50um_25hz', 'EWS_25um_25hz', 'circ_1000hz', 'sqr_1000hz', 'P_1000hz', 'T_1000hz', 'E_1000hz', 'X_1000hz'};
testToCheck = 7; % test number required to correspond to the sorted cell data
numCells = max(size(cacahuete.listOfCells_spk)); %number of cells

for cells = 1:numCells
    index = find(cacahuete.listOfCells_comeFrom{cells} == testToCheck); % separate spikes corresponding to the correct test
    index = index(1:2:end);
    spikesToConsider{cells} = double(cacahuete.listOfCells_spk{cells}(index)*SAMPLERATE);   
    channelsOfSpike(cells) = cacahuete.listOfCells_elec{cells}(1);
end


%% load spike times
load([pathMain '/cellWiseAxonShapes.mat']);
%% load spike times
% startOfStimulation = rawData(1);
% endOfStimulation = frameTimes(end);


% for each cell get the raw signal of all the electrodes when this cell
% spikes

disp('loading spike data..');
shapeSize = 101;

%The variable targetShapes contains mean signal shape of one cell w.r.t all other cells

hiFreq = 1000;
[b,a] = butter(2, hiFreq/(sampling_rate/2), 'high');


% targetShape = zeros(numCells, numCells, shapeSize);



parfor eachCell = 1:numCells %for each spike of the reference cell
    refData = double(m.Data(channelsOfSpike(eachCell):NUM_CHANNELS:NUM_CHANNELS*SAMPLERATE*Time));
    refData = (refData - h_r.adc_zero)*h_r.conversion_factor;
    refData =    filter(b,a,refData);
    
    SpikeTimes = find(refData(1:end-1)>-15 & refData(2:end)<=-15);
    spikeTimesRef = SpikeTimes; %round(spikesToConsider{eachCell}); 

    spikeTimesRef = spikeTimesRef(spikeTimesRef > 50);
    
    disp(eachCell)
    tic
    for eachOtherCell =  1:numCells
        targetData = double(m.Data(channelsOfSpike(eachOtherCell):NUM_CHANNELS:NUM_CHANNELS*h_r.npts));
        targetData = (targetData - h_r.adc_zero)*h_r.conversion_factor;
        targetData =    filter(b,a,targetData);
        
        temp = zeros(numel(spikeTimesRef)-3, 101);
        for eachSpikeOfCell = 1:numel(spikeTimesRef)-3 %spikesToConsider{eachCell})
%             temp = squeeze(targetShape(eachCell, eachOtherCell, :));
            
            %get 100 points => 5 ms of raw data around the spike
%             temp = temp + targetData( floor(spikesToConsider{eachCell}(eachSpikeOfCell)) - 50: floor(spikesToConsider{eachCell}(eachSpikeOfCell)) + 50);
            temp(eachSpikeOfCell,:) = targetData( floor(spikeTimesRef(eachSpikeOfCell)) - 50: floor(spikeTimesRef(eachSpikeOfCell)) + 50);
           
%             targetShape(eachCell, eachOtherCell, :) = temp;
%             [eachCell eachOtherCell eachSpikeOfCell]
        end
            %mean the spike shape over all spikes
            targetShape(eachCell, eachOtherCell, :) = mean(temp); %targetShape(eachCell, eachOtherCell, :) / (numel(spikeTimesRef)-3 );
    end
    toc
end

%% Compute latencies
disp('computing latencies...');

%latency is given as the difference between the minimum point of the signal


latency = zeros(numCells, numCells);
for eachCells =  1:numCells
    
        refShape = squeeze(targetShape(eachCells, eachCells, :));
        refShape = refShape - mean(refShape);
        
        for eachOtherCells = 1:numCells
            spikeShape = squeeze(targetShape(eachCells, eachOtherCells, :));
            spikeShape = spikeShape - mean(spikeShape);
            

                latency(eachCells, eachOtherCells) = (find(spikeShape == min(spikeShape), 1) - find(refShape == min(refShape), 1) )/20 ; %convert index to milli-seconds [sampling freq = 20000 hz]

        end
end

%% Compute speeds
disp('computing axonal speeds');
speedComputes = [];


% opto
MEA_MAP1 = MEA_MAP';
MEA_MAP1 = MEA_MAP1(:);

% multipix
% MEA_MAP1 = MEA_MAP;
n=0;

% gcf = figure(eachCells); 
%         set(gcf, 'Position', [0 0 1000 1000]);
%         imagesc(image1); axis off; axis square;
%         colormap(gcf, hot);
%         hold on; 
% 
%         drawnow;
        
for eachCells = 1:numCells
    
        refShape = squeeze(targetShape(eachCells, eachCells, :));
        refShape = refShape - mean(refShape); % normalize the spike shape so that the signal is mostly zero accept the spike area
        c=0;
        lats = [];

        for eachOtherCells = 1:numCells %channelsToConsider
            spikeShape = squeeze(targetShape(eachCells, eachOtherCells, :));
            spikeShape = spikeShape - mean(spikeShape); % normalize the spike shape so that the signal is mostly zero accept the spike area
            
            % for a pair of cells/electrodes compute the correlation
            % between their shapes (crosscorr). If they are similar only
            % then consider them on an axon
            if(min(spikeShape)/min(refShape) > 0.01 && max((crosscorr(spikeShape, refShape))) > 0.9 )
                
                %if the cells/electrodes are on same axons then 
                
                % 1. Compute the latency between the ref spike and the
                % target electrode spike
                
                
                c=c+1;
                lats(c) = latency(eachCells, eachOtherCells);
                cellsLats(c) = eachOtherCells;
            
                
                % 2. Find where this electrode is located. 
                % This is used to compute the distance and to plot the
                % axons 
                indexnChannel = cellsLats(c);
                cChannel = channelsOfSpike(indexnChannel(1));
                cChannel = find(cChannel == MEA_MAP1);            
                cChannelPos = electrode(cChannel).coord_DMD;
            
                cellLocations(c,:) = cChannelPos;
                
                % The variable 'c' is to count the number of electrodes/cells found on an axon.
                % This is used to have another confidence level to ensure
                % we are considering an axon
            end
        end
        
%         speedComputes(eachCells,1 ) = NaN;

        % if the variable 'c' is more than 6 i.e. there are atleast 6
        % cells/electrodes on an axon, we will compute the speeds to avoid
        % spurious computations 
        
        
        % The variable speedComputes contains all the information about the
        % axons such as:
        % 1. The total latency between the fastest and slowest spike
        % incidence (Lat) .. in msec
        % 2. Total distance between the farthest and closest electode
        % (Distance) .. in  um 
        % 3. The axon speed: given as Distance/Lat (we multiply 0.001 to
        % this to convert from um/msec to m/sec
        % 4. The individual values of the smallest, largest latencies and
        % electrodes and the position of these electrodes
        % 
        if(c > 6)
            n=n+1;
            
           
            
            speedComputes(n, 1) = max(lats) - min(lats); % total latency difference
            
            indexMinChannel = find(lats == min(lats));
            indexMinChannel = cellsLats(indexMinChannel); % find the low latency electrode
            
            minChannel = channelsOfSpike(indexMinChannel(1));
            minChannel = find(minChannel == MEA_MAP1);            
            minChannelPos = electrode(minChannel).coord_DMD; % find the low latency electrode location
            
            indexMaxChannel = find(lats == max(lats));
            indexMaxChannel = cellsLats(indexMaxChannel); % find the highest latency electrode
            
            maxChannel = channelsOfSpike(indexMaxChannel(1));
            maxChannel = find(maxChannel == MEA_MAP1);
            maxChannelPos = electrode(maxChannel).coord_DMD; % find the highest latency electrode location
            
            speedComputes(n, 5) = minChannel;
            speedComputes(n,6) = maxChannel;
            
            speedComputes(n, 2) = norm(minChannelPos - maxChannelPos)*(2.23); % Total distance between the highest and lowest latency electrode = Axon Length
            
%             netLatencies = netLatencies(netLatencies ~=0);
%             lengthAxon = lengthAxon(lengthAxon ~=0);
            
            speedComputes(n, 3) = 0.001*speedComputes(n, 2)/speedComputes(n, 1); % Axon speed. % distance is in um and time in milli-second - so multiply by 1000 to convert to meter/sec
            speedComputes(n, 4) = eachCells;
 
            speedComputes(n, 7) = min(lats);
            speedComputes(n, 8) = max(lats);
            
            speedComputes(n, 9) = minChannelPos(1);
            speedComputes(n, 10) = minChannelPos(2);
            
            speedComputes(n, 11) = maxChannelPos(1);
            speedComputes(n, 12) = maxChannelPos(2);
            
             if(speedComputes(n,4)~=128 && speedComputes(n,3)<=3.5 && speedComputes(n,4)~=127)
                [sortedLats, indexSorted] = sort(lats, 'ascend');
                speedComputesAllLats{n} = sortedLats;
                tempCellLocations = cellLocations(indexSorted,:);

                tempDistCompute = tempCellLocations - repmat(tempCellLocations(1,:), size(tempCellLocations,1), 1); %[502.8646 375.1736]
                tempDistCompute = sqrt(tempDistCompute(:, 1).^2 + tempDistCompute(:, 2).^2);

                speedComputesAllLocs{n} = tempDistCompute;
             end
            
%             for tempI_ = 1:c
%                 indexTempChannel = cellsLats(tempI_); %find(latency(eachCells,:) == lats(tempI_));
%                 tempChannel = channelsOfSpike(indexTempChannel(1));
%                 
%                 k = find(tempChannel == MEA_MAP1);
%                 hold on;
%                 pointChannel = [electrodeLocations{k}(2), electrodeLocations{k}(1)];
% 
%                 scatter(pointChannel(1), pointChannel(2), 30, 'b', 'Filled');                    
%                 text(pointChannel(1), pointChannel(2), [num2str(tempChannel)], 'Fontsize', 10, 'FontWeight', 'Bold', 'color', 'g');
%             end
        
        end
        

end
% hist(speedComputes, 0.1:0.1:1.6);

% As a precaution, remove the axons if the electrode was found on channel
% 127 or 128 
speedComputes1 = speedComputes(speedComputes(:,4) ~= 128,:);
% speedComputes1 = speedComputes1(speedComputes1(:,3) <=2.5,:);
speedComputes1 = speedComputes1(speedComputes1(:,4) ~= 127,:);
speedComputes = speedComputes1;

% convert the matrix into a data set with column names as indicated.
speedComputes = dataset({speedComputes 'LatDiff' 'ElectrodeDist' 'AxonSpeed' 'CellNum' 'minChLoc' 'maxChanLoc' 'minLat' 'maxLat' 'minChPos1' 'minChPos2' 'maxchPos1' 'maxChPos2'});


%%

% MEA_MAP1 = MEA_MAP(:, end:-1:1);
MEA_MAP1 = MEA_MAP';
MEA_MAP1 = MEA_MAP1(:);


adjustX = 0;
adjustY = -0;
rot = 0.0;

for e=1:256
    
    P(1) = electrode(e).coord_TIFF(1);% + adjustX - electrode(16).coord_TIFF(1);
    P(2) = electrode(e).coord_TIFF(2);% + adjustY - electrode(16).coord_TIFF(2);
    
%     P = P*[cos(rot) sin(rot); -sin(rot) cos(rot)] + [electrode(16).coord_TIFF(1) electrode(16).coord_TIFF(2)];
%     hold on; scatter(P(1) , P(2), 80, 'y', 'Filled');
%     text(P(1), P(2), num2str(MEA_MAP(e)), 'Color', cols(100,:), 'FontWeight', 'Bold', 'FontAngle', 'italic');
    electrodeLocations{e} = P;
end


%%
plotFigures = true;
if(plotFigures)

% Get the axon data and plot the axons and the mean spike shapes on the
% florescence image

%     close all;

    [image, ~] = imread([pathMain '/test0f.tif']);

    
    image1 = double(squeeze(image(:,:,1)));
    % image1 = image1(:, end:-1:1);
    image1 = image1/max(max(image1));
    % image1 = imrotate(image1, 10, 'crop');
    % image1(:,:, 1~) = squeeze(image(:,:,1));

    mkdir([pathDir '/AxonFinder_jpgs/']);

    gcf = figure(1); 
    set(gcf, 'Position', [0 0 1000 1000]);
    
    for eachCells = 1:100 %numCells
            eachCells
    %     if(channelWiseSTA.numSpikes(channel) > 1000)
            
            flag=0;
            c=0;
            refShape = squeeze(targetShape(eachCells, eachCells, :));
            refShape = refShape - mean(refShape);
            [~, refTime] = max(refShape);
            
            delayes = [];
            shapeToPlot = [];
            channelToPlot = [];
            isSameChannel = [];

            for eachOtherCells = 1:numCells
                eachOtherChannel = channelsOfSpike(eachOtherCells);

                spikeShape = squeeze(targetShape(eachCells, eachOtherCells, :));
                spikeShape = spikeShape - mean(spikeShape);
                if(max(spikeShape)/max(refShape) > 0.19 && eachOtherChannel~=15 )
                    flag=1;
                    c=c+1;
                    [~, targetTime] = max(targetShape(eachCells, eachOtherCells, :));
                    
                    delayes(c) = refTime - targetTime;
                    shapeToPlot(c, :) = spikeShape;
                    channelToPlot(c) = eachOtherChannel;
                    cellToCheck(c) = eachOtherCells;
                    
                    isSameChannel(c) = (eachOtherCells == eachCells);
                end
            end
            
            if(c>6)
                [eachCells c]
                [delayes_, indxSorted] = sort(delayes);                
                
                                
                delayes_ = delayes_ - min(delayes_) + 1;
%                 delayes_(12) = [];
%                 delayes_(25) = [];
                
                
                delayes = delayes_;
                
                cols = hcparula(max(delayes)+5);
               
                imagesc(image1); axis off; axis square;
                colormap(gcf, gray);            
                hold on; 
               
%                             [1  2  3  4  5  6  7  8  9  10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27];
%                 c_display = [26 27 24 25 23 22 21 20 19 17 18  6 16 14  5 11 15 12 13  3 10  4  2  9  8  1  7];
%                 c_display =  [27 26 25 24 23 22 21 20 19 18 17 16 15 14 13 12 11 10 9 8 7 6 5 4 3 2 1];
                
                c_display = 1:numel(delayes);
                for c_ = 1:numel(delayes) %eachOtherCells = 1:numCells[1:11 13:numel(delayes)-1] %[1:11 12 13:numel(delayes)-2 numel(delayes) ] 
                    
                    k = find(channelToPlot(indxSorted(c_)) == MEA_MAP1);
                    pointChannel = [electrodeLocations{k}(2), electrodeLocations{k}(1)];                        
                    toPlot = shapeToPlot(indxSorted(c_),:);
                    
                    
                    if(isSameChannel(indxSorted(c_)))
                        %                     scatter(pointChannel(1), pointChannel(2), 90, 'r', 'Filled');
                        ii = -50:50;
%                         scatter(pointChannel(1)-ii, pointChannel(2)-20*toPlot(ii+51), 30, cols(delayes(c_)+5,:), 'Filled');
                        scatter(pointChannel(1), pointChannel(2), 600, cols(delayes(c_)+5,:), 'Filled');
%                         scatter(pointChannel(1), pointChannel(2), 30, 'b', 'Filled');
%                         for i = 0%-50:50                        
%                             scatter(pointChannel(1)-i, pointChannel(2), 30, 'b', 'Filled');
%                             %                         scatter(pointChannel(1)-i, pointChannel(2)-50*toPlot(i+51), 30, 'b', 'Filled');
%                         end
                    else
                        ii = -50:50;
%                         scatter(pointChannel(1)-ii, pointChannel(2)-20*toPlot(ii+51), 30, cols(delayes(c_)+5,:), 'Filled');
                        scatter(pointChannel(1), pointChannel(2), 600, cols(delayes(c_)+5,:), 'Filled');
%                         for i = -50:50
%                             %                             scatter(pointChannel(1)-i, pointChannel(2), 30, 'y', 'Filled');                        
%                             scatter(pointChannel(1)-i, pointChannel(2)-20*toPlot(i+51), 30, cols(delayes(c),:), 'Filled');
%                         end
                    end
%                     c_display = c_display  +1;
                    if(c_display(c_) < 10)
                        text(pointChannel(1)-5, pointChannel(2), num2str(c_display(c_) ), 'Color', 'm', 'FontSize', 10, 'FontWeight', 'bold');
%                         text(pointChannel(1)-5, pointChannel(2)-15, num2str(c_), 'Color', 'k', 'FontSize', 10, 'FontWeight', 'bold');
                    else
                        text(pointChannel(1)-8, pointChannel(2), num2str(c_display(c_) ), 'Color', 'm', 'FontSize', 10, 'FontWeight', 'bold');
%                         text(pointChannel(1)-5, pointChannel(2)-15, num2str(c_), 'Color', 'k', 'FontSize', 10, 'FontWeight', 'bold');
                    end
                    %                 plot(channelWiseSTA.rawSpikes(:, channel, eachOtherChannel)/channelWiseSTA.numSpikes(channel));
                    %                 xlim([0 101]); ylim([-0.02 0.02]);
                    drawnow;
                end
                
%                 saveas(gcf, [pathDir '/AxonFinder_jpgs/' date_ '-10-2016-channel-stim-' num2str(spotSize) 'um_' num2str(freq) 'Hz_axons_' num2str(eachCells) '_shapes'], 'png')
%                 saveas(gcf, [pathDir '/AxonFinder_jpgs/' date_ '-10-2016-channel-stim-' num2str(spotSize) 'um_' num2str(freq) 'Hz_axons_' num2str(eachCells) '_shapes'], 'fig')
                
%                 saveas(gcf, [pathDir '/AxonFinder_jpgs/' date_ '-10-2016-channel-stim-' num2str(spotSize) 'um_' num2str(freq) 'Hz_axons_' num2str(eachCells) '_dots_sized'], 'png')
%                 saveas(gcf, [pathDir '/AxonFinder_jpgs/' date_ '-10-2016-channel-stim-' num2str(spotSize) 'um_' num2str(freq) 'Hz_axons_' num2str(eachCells) '_dots_sized'], 'fig')

                saveas(gcf, [pathDir '/AxonFinder_jpgs/' date_ '-10-2016-channel-stim-' num2str(spotSize) 'um_' num2str(freq) 'Hz_axons_' num2str(eachCells) '_dots_new'], 'png');
                saveas(gcf, [pathDir '/AxonFinder_jpgs/' date_ '-10-2016-channel-stim-' num2str(spotSize) 'um_' num2str(freq) 'Hz_axons_' num2str(eachCells) '_dots_new'], 'fig');

%                 saveas(gcf, [pathDir '/AxonFinder_jpgs/' date_ '-10-2016-channel-stim-' num2str(spotSize) 'um_' num2str(freq) 'Hz_axons_' num2str(eachCells) '_dots'], 'svg')
            end
clf(gcf);            

    %     end
    end
    close all;
end

%% Analyze raw data shapes for supp figure


% startOfStimulation = rawData(1);
% endOfStimulation = frameTimes(end);


% for each cell get the raw signal of all the electrodes when this cell
% spikes

disp('loading spike data..');
shapeSize = 101;

%The variable targetShapes contains mean signal shape of one cell w.r.t all other cells

targetShape = zeros(numCells, numCells, shapeSize);

close all;

for eachCell = 1:100%numCells %for each spike of the reference cell
    refData_ = double(m.Data(channelsOfSpike(eachCell):NUM_CHANNELS:NUM_CHANNELS*SAMPLERATE*Time));
  %%
    refData = (refData_ - h_r.adc_zero)*h_r.conversion_factor;
    hiFreq = 1000;
    [b,a] = butter(2, hiFreq/(sampling_rate/2), 'high');
    
    refData = filter(b,a,refData);  
    SpikeTimes = find(refData(1:end-1)>-15 & refData(2:end)<=-15);
    spikeTimesRef = SpikeTimes; %round(spikesToConsider{eachCell}); 

  %%
    close all;
    length_ = 50;

    c__=1;
    
    gcf = figure(1);
    temp = [];
    for n = 1:10 %numel(spikeTimesRef) %[7 8 9 10 23 40 43 56 57 58 66 70] 
        
        
%         valueSpike(n) = refData(spikeTimesRef(n));
%         if(refData(spikeTimesRef(n)) < -10)
            
            plot(c__:c__+numel(refData(spikeTimesRef(n)-length_ :spikeTimesRef(n)+length_ ))-1, refData(spikeTimesRef(n)-length_ :spikeTimesRef(n)+length_ ), 'k', 'LineWidth', 2);
            
            temp = [temp refData(spikeTimesRef(n)-length_ :spikeTimesRef(n)+length_ )];
            hold on; 
            plot(c__+numel(refData(spikeTimesRef(n)-length_ :spikeTimesRef(n)+length_ ))/2, 10, '.r', 'MarkerSize', 16);
%             text(c__+numel(refData(spikeTimesRef(n)-length_ :spikeTimesRef(n)+length_ ))/2, 15, num2str(n), 'Color', 'k', 'FontSize', 6);
            
            c__ = c__ + numel(refData(spikeTimesRef(n)-length_ :spikeTimesRef(n)+length_ ));            
%         end
    end
            plot(c__ + 50:c__ + 50 + numel(mean(temp'))-1, 2*mean(temp'), 'Color', 'k', 'LineWidth', 2);

    
%     saveas(gcf, [pathDir '/AxonFinder_jpgs/' date_ '-10-2016-channel-stim-' num2str(spotSize) 'um_' num2str(freq) 'Hz_axons_' num2str(eachCells) '_raw_signals'], 'png')
%     saveas(gcf, [pathDir '/AxonFinder_jpgs/' date_ '-10-2016-channel-stim-' num2str(spotSize) 'um_' num2str(freq) 'Hz_axons_' num2str(eachCells) '_raw_signals'], 'fig')      
    
%     clf(gcf);
    
    c1 = 2;
    colmaps_ = hcparula(numel(delayes)+3);
%     colmaps_ = colmaps_(1:numel(delayes), :);
    for n__ = 1:size(colmaps_,1)
        colmaps(size(colmaps_,1)-n__+1,:) = colmaps_(n__, :);
    end
%     colmaps = colmaps_;

    
    for eachOtherCell =  1:numCells %[49 53 76 54 34 14 77 ] %[70 79 35 78 51 70]
        % [60 38 61 82 80 59 79 16 35 55 56 78 51 2 50 72 7 24 27 48 4 70]
        % [2 4 7 16 24 27 35 38 48 50 51 55 56 59 60 61 70 72 78 79 80 82]%1:numCells
        c1 = c1+1;
        targetData_ = double(m.Data(channelsOfSpike(eachOtherCell):NUM_CHANNELS:NUM_CHANNELS*h_r.npts));        
        targetData = (targetData_ - h_r.adc_zero)*h_r.conversion_factor;
        targetData = filter(b,a,targetData);
        
        temp = [];
        
        c__=1;
        for n = 1:10%numel(spikeTimesRef) %[7 8 9 10 23 40 43 56 57 58 66 70]    %1:numel(spikeTimesRef) %    
%             if(refData(spikeTimesRef(n)) < -10)
                
                temp = [temp targetData(spikeTimesRef(n)-length_ :spikeTimesRef(n)+length_ )];

%                 plot(c__:c__+numel(targetData(spikeTimesRef(n)-length_ :spikeTimesRef(n)+length_ ))-1,  30*c1 + targetData(spikeTimesRef(n)-length_ :spikeTimesRef(n)+length_ ), 'Color', colmaps(c1,:), 'LineWidth', 2);
                plot(c__:c__+numel(targetData(spikeTimesRef(n)-length_ :spikeTimesRef(n)+length_ ))-1,  30*c1 + targetData(spikeTimesRef(n)-length_ :spikeTimesRef(n)+length_ ), 'Color', 'k', 'LineWidth', 2);
                
                hold on; 
%                 plot(c__+numel(targetData(spikeTimesRef(n)-200:spikeTimesRef(n)+200))/2, 10, '.r', 'MarkerSize', 16);
                
                c__ = c__ + numel(targetData(spikeTimesRef(n)-length_ :spikeTimesRef(n)+length_ ));            
%             end
        end        
        drawnow;
        
        plot(c__ + 50:c__ + 50 + numel(mean(temp'))-1, 30*c1 + 2*mean(temp'), 'Color', 'k', 'LineWidth', 2);
        %numel(squeeze(targetShape(eachCell, eachOtherCell, :)))-1
        
        
%         saveas(gcf, [pathDir '/AxonFinder_jpgs/' date_ '-10-2016-channel-stim-' num2str(spotSize) 'um_' num2str(freq) 'Hz_axons_' num2str(eachCells) '_raw_signals_' num2str(eachOtherCell) ], 'png')
%         saveas(gcf, [pathDir '/AxonFinder_jpgs/' date_ '-10-2016-channel-stim-' num2str(spotSize) 'um_' num2str(freq) 'Hz_axons_' num2str(eachCells) '_raw_signals_'  num2str(eachOtherCell)], 'fig')      
%         clf(gcf);
    end    
        c__ = 1;
        for n = 1:10% numel(spikeTimesRef) %[7 8 9 10 23 40 43 56 57 58 66 70] %1:numel(spikeTimesRef)        
%             if(refData(spikeTimesRef(n)) < -10)
                hold on;
                patch([c__+(length_)-30, c__+(length_)-30, c__+(length_)+30, c__+(length_)+30], [-30 (c1+1)*30 (c1+1)*30 -30], [0.5 0.5 0.5], 'FaceAlpha', 0.2);
                c__ = c__ + numel(targetData(spikeTimesRef(n)-length_ :spikeTimesRef(n)+length_ ));           
%             end
        end
    
%     saveas(gcf, [pathDir '/AxonFinder_jpgs/' date_ '-10-2016-channel-stim-' num2str(spotSize) 'um_' num2str(freq) 'Hz_axons_' num2str(eachCells) '_raw_signals_all' ], 'png')
%     saveas(gcf, [pathDir '/AxonFinder_jpgs/' date_ '-10-2016-channel-stim-' num2str(spotSize) 'um_' num2str(freq) 'Hz_axons_' num2str(eachCells) '_raw_signals_all'], 'fig')      
end
















