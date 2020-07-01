
%%% THIS CODE IS INCOMPLETE AND SHALL NOT BE USED YET


%% initialize params
clear

NUM_CHANNELS = 256;
SAMPLERATE = 20000;
totalTimeStim = 20;
timeStim = 10;
reps = 25;

spotSize = 25;

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

%%
sampling_rate = 20000;

%% load vec and frametime file

date_ = '11';
pathMain = ['../../../DATA/' date_ '_10_2016_opto_primate/'];

load([pathMain '/' 'cacahuete_day_20161011_EWS__.mat']);

[file, pathDir] = uigetfile([pathMain '/*.raw'],'Select any channel file:' , 'MultiSelect', 'on');
indexEnd = strfind(file, '_C_');

%%

f = size(file, 1);

THRESHOLD = 50;

h_r = mcd_header([pathDir file]);

m = memmapfile([pathDir file], 'Offset', h_r.data_start, 'Format', 'uint16');


% h_r = mcd_header([pathDir Filename{f} '.raw']);
% 
% m = memmapfile([pathDir Filename{f} '.raw'], 'Offset', h_r.data_start, 'Format', 'uint16');

Time = h_r.npts/SAMPLERATE;


[b,a] = butter(2, 300/(SAMPLERATE /2), 'high');

%%

testToCheck = 9;
numCells = max(size(cacahuete.listOfCells_spk));

for cells = 1:numCells
    index = find(cacahuete.listOfCells_comeFrom{cells} == testToCheck);
    index = index(1:2:end);
    spikesToConsider{cells} = double(cacahuete.listOfCells_spk{cells}(index)*SAMPLERATE);   
    channelsOfSpike(cells) = cacahuete.listOfCells_elec{cells}(1);
end

%%
referenceShape = zeros(numCells, 101);

for eachCell = 1:numCells
    data = double(m.Data(channelsOfSpike(eachCell):NUM_CHANNELS:NUM_CHANNELS*SAMPLERATE*Time));
    data = (data - h_r.adc_zero)*h_r.conversion_factor;
        
    for eachSpikeOfCell = 1:numel(spikesToConsider{eachCell})
%         if(~isinteger(spikesToConsider{eachCell}(eachSpikeOfCell)) )
%             disp('why');
%             spikesToConsider{eachCell}(eachSpikeOfCell)
%             continue;
%         end
        referenceShape(eachCell, :) = referenceShape(eachCell, :) + data( floor(spikesToConsider{eachCell}(eachSpikeOfCell)) - 50: floor(spikesToConsider{eachCell}(eachSpikeOfCell)) + 50)';
    end
    referenceShape(eachCell, :) = referenceShape(eachCell, :) / numel(spikesToConsider{eachCell});
end

%%

% for channel = ChannelsToConsider
%     
%     dataStart = channel;
% 
%     data = double(m.Data(dataStart:NUM_CHANNELS:NUM_CHANNELS*sampling_rate*Time));
%     data = (data - h_r.adc_zero)*h_r.conversion_factor;
% 
% 	rawData{channel} = data;   
% end



%% load spike times

% startOfStimulation = rawData(1);
% endOfStimulation = frameTimes(end);

shapeSize = 101;

targetShape = zeros(numCells, numCells, shapeSize);

for eachCell = 1:numCells
%     refData = double(m.Data(channelsOfSpike(eachCell):NUM_CHANNELS:NUM_CHANNELS*SAMPLERATE*Time));
%     refData = (redData - h_r.adc_zero)*h_r.conversion_factor;
        
    for eachOtherCell = 1:numCells
        targetData = double(m.Data(channelsOfSpike(eachOtherCell):NUM_CHANNELS:NUM_CHANNELS*SAMPLERATE*Time));
        targetData = (targetData - h_r.adc_zero)*h_r.conversion_factor;
        for eachSpikeOfCell = 1:numel(spikesToConsider{eachCell})
            temp = squeeze(targetShape(eachCell, eachOtherCell, :));
            temp = temp + targetData( floor(spikesToConsider{eachCell}(eachSpikeOfCell)) - 50: floor(spikesToConsider{eachCell}(eachSpikeOfCell)) + 50);
           
            targetShape(eachCell, eachOtherCell, :) = temp;
%             [eachCell eachOtherCell eachSpikeOfCell]
        end
            targetShape(eachCell, eachOtherCell, :) = targetShape(eachCell, eachOtherCell, :) / numel(spikesToConsider{eachCell});
    end

end

%%
% 
% channelWiseSTA.rawSpikes = zeros(101, 256, 256);
% channelWiseSTA.numSpikes(channel) = 0;
% 
%         
% for cells = 1:numCells %ChannelsToConsider
% %         spikedata = load([pathDir '/ST_' file(8:indexEnd+2) num2str(channel) '.mat']);   
% %         spikeTimes = spikedata.SpikeTime;
%         channel
%         
%         spikeTimesChannel = spikeData{channel}.spikeTimes;
% %         rawChannel = spikeData{channel}.raw;
% 
%         numel(spikeTimesChannel)
% 
%         for eachSpike = 1:numel(spikeTimesChannel)
%             for eachOtherChannel = ChannelsToConsider
% 
% %                load([pathDir '/' file(1:indexEnd+2) num2str(eachOtherChannel) '.mat']);   
%                 raw = spikeData{eachOtherChannel}.raw;
% 
% 	
%                 if( (spikeTimesChannel(eachSpike)-50) > 0 && (spikeTimesChannel(eachSpike)+50) < numel(raw))
%                     channelWiseSTA.numSpikes(channel) = channelWiseSTA.numSpikes(channel)+1;
%                     channelWiseSTA.rawSpikes(:, channel, eachOtherChannel) = channelWiseSTA.rawSpikes(:, channel, eachOtherChannel) + raw(spikeTimesChannel(eachSpike)-50:spikeTimesChannel(eachSpike)+50); 
%                 end
%             
%             end
%         end
% end

%%

% MEA_MAP1 = MEA_MAP(:, end:-1:1);
MEA_MAP1 = MEA_MAP';
MEA_MAP1 = MEA_MAP1(:);

cols = summer(20);

adjustX = 22;
adjustY = -60;
rot = 0.0;

for e=1:256
    
    P(1) = electrode(e).coord_TIFF(1) + adjustX - electrode(16).coord_TIFF(1);
    P(2) = electrode(e).coord_TIFF(2) + adjustY - electrode(16).coord_TIFF(2);
    
    P = P*[cos(rot) sin(rot); -sin(rot) cos(rot)] + [electrode(16).coord_TIFF(1) electrode(16).coord_TIFF(2)];
%     hold on; scatter(P(1) , P(2), 80, 'y', 'Filled');
%     text(P(1), P(2), num2str(MEA_MAP(e)), 'Color', cols(100,:), 'FontWeight', 'Bold', 'FontAngle', 'italic');
    electrodeLocations{e} = P;
end


%%
close all;

% [image, ~] = imread([pathMain '/2016_04_20_BA901h_OD1_tdt_4x_4s_13.6gain.tif']);
[image, ~] = imread([pathMain '/test0f.tif']);

image1(:,:, 1) = squeeze(image(:,:,1));


for channel = ChannelsToConsider
    
%     if(channelWiseSTA.numSpikes(channel) > 1000)
        gcf = figure(channel); 
        set(gcf, 'Position', [0 0 1000 1000]);
        imagesc(image1/200); axis off; axis square;
        colormap(gcf, gray);
        hold on; 

        flag=0;
        c=0;
        for eachOtherChannel = ChannelsToConsider
            
            if(max(abs(channelWiseSTA.rawSpikes(:, channel, eachOtherChannel)/channelWiseSTA.numSpikes(channel)) - mean(abs(channelWiseSTA.rawSpikes(:, channel, eachOtherChannel)/channelWiseSTA.numSpikes(channel))) ) > 0.01)
                flag =1;
                
                c=c+1;
%                 figure(channel);
%                 set(gcf, 'Position', [0, 0, 1000, 1000]);
%                 drawnow;
%                 subplot(16, 16, MEA_MAP1(eachOtherChannel));

%                 pointElectrode = [electrodeLocations{channel}(1), electrodeLocations{channel}(2)];

                k = find(eachOtherChannel == MEA_MAP1);
                pointChannel = [electrodeLocations{k}(1), electrodeLocations{k}(2)];

                toPlot = channelWiseSTA.rawSpikes(:, channel, eachOtherChannel)/channelWiseSTA.numSpikes(channel) - mean(channelWiseSTA.rawSpikes(:, channel, eachOtherChannel)/channelWiseSTA.numSpikes(channel));
                if(eachOtherChannel == channel)
%                     scatter(pointChannel(1), pointChannel(2), 90, 'r', 'Filled');

                    for i = -50:50                        
                        scatter(pointChannel(1)-i, pointChannel(2)-5000*toPlot(i+51), 30, 'b', 'Filled');
                    end
                else
                    for i = -50:50
                        scatter(pointChannel(1)-i, pointChannel(2)-5000*toPlot(i+51), 30, cols(c,:), 'Filled');
                    end
                end
%                 plot(channelWiseSTA.rawSpikes(:, channel, eachOtherChannel)/channelWiseSTA.numSpikes(channel));
%                 xlim([0 101]); ylim([-0.02 0.02]);
            end
        end
        if(c>6)
                    drawnow;
            saveas(gcf, [pathDir '/jpgs/' date_ '-10-2016-channel-stim-' num2str(spotSize) 'um_' num2str(freq) 'Hz_axons_' num2str(channel) ], 'png')
            saveas(gcf, [pathDir '/jpgs/' date_ '-10-2016-channel-stim-' num2str(spotSize) 'um_' num2str(freq) 'Hz_axons_' num2str(channel) ], 'fig')
        end
        close all;
        
            
%     end
end

