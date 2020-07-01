%% initialize params
% clear


sampleRate = 20000;
totalTimeStim = 20;
timeStim = 10;
reps = 25;

spotSize = 100;

freq = 100;%floor(1000/timeStim); %floor(1000/totalTimeStim);

considerWindowMs = 500; %ms Window to consider post stim
considerWindowInd = considerWindowMs*sampleRate/1000; %num_points Window to consider post stim

time = -considerWindowMs:considerWindowMs; %ms Window to consider post stim


ChannelsToConsider = [2:14 16:126 129:254];


rasterData = cell(256, 256);


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

% MEA_MAP1 = MEA_MAP(end:-1:1, :);
MEA_MAP1 = MEA_MAP(:, end:-1:1);
% MEA_MAP1 = MEA_MAP1';
MEA_MAP1 = MEA_MAP1(:);
        
%% load vec and frametime file

vec = load(['C:\Users\himanshu-idv\POST_DOC\DATA\MEA_DATA\BINs_VECs\tps_' num2str(timeStim) '_reps_' num2str(reps) '_periode_' num2str(totalTimeStim) '.vec']);
vec = vec(2:end,2);

pathMain = 'C:\Users\himanshu-idv\POST_DOC\DATA\MEA_DATA\Primates\20_04_2016\';
pathDir = ['C:\Users\himanshu-idv\POST_DOC\DATA\MEA_DATA\Primates\20_04_2016\SpikeTimes\2016_04_20_BA901h_OD1_gsmmix_DMD600i100_EWS_' num2str(spotSize) 'um_' num2str(timeStim) 'ms_' num2str(freq) 'hz/'];

frameTimes = load([pathDir '/ST_2016_04_20_BA901h_OD1_gsmmix_DMD600i100_EWS_' num2str(spotSize) 'um_' num2str(timeStim) 'ms_' num2str(freq) 'hz_C_127.mat']);   
frameTimes = frameTimes.SpikeTime;


load('C:\Users\himanshu-idv\POST_DOC\DATA\MEA_DATA\Primates\20_04_2016\test6_20-Apr-2016.mat');      


%% load spike times
for eachElectrode = 1:256
    eachElectrode
    
    indStim = find(vec == eachElectrode);    
    numTrials = numel(indStim)/(timeStim);
    
    indStimTemp =  find( mod((indStim-1), totalTimeStim) == 0);
    
%     if(sum(indStim(indStimTemp) > numel(frameTimes)))
%         continue;
%     end
            
    tempIndices = indStim(indStimTemp); %indStim( indStim(indStimTemp) < numel(frameTimes));
    startOfStims = frameTimes(tempIndices);
    
    indStimTemp =  find(mod((indStim), totalTimeStim+timeStim) == 0);
    tempIndices = indStim( indStim(indStimTemp) < numel(frameTimes));
    endOfStims = frameTimes(tempIndices);
    
    maxTrials(eachElectrode) = numel(startOfStims);

    for channel = ChannelsToConsider
    spikeTimes = load([pathDir '/ST_2016_04_20_BA901h_OD1_gsmmix_DMD600i100_EWS_' num2str(spotSize) 'um_' num2str(timeStim) 'ms_' num2str(freq) 'hz_C_' num2str(channel) '.mat']);   
    spikeTimes = spikeTimes.SpikeTime;

    rasterData{eachElectrode, channel} = cell(maxTrials(eachElectrode),1);
        for trial = 1:numel(startOfStims)
                tempSpikes = spikeTimes( spikeTimes > (startOfStims(trial) - considerWindowInd) & spikeTimes < (startOfStims(trial) + considerWindowInd) ) - startOfStims(trial);
%                 if(~isempty(tempSpikes)) 
%                     [trial channel]
%                 end
                rasterData{eachElectrode, channel}{trial} = tempSpikes; %spikeTimes( spikeTimes > startOfStims(trial) & spikeTimes < startOfStims(trial) + considerWindowInd) - startOfStims(trial);
        end
    end
end

%% compute spike density response
close all;
sdf = cell(256, 256);

for e = 1:256
    e
    for c = ChannelsToConsider
        for t = 1:maxTrials(eachElectrode)
            temp = rasterData{e,c}{t};
            
            if(~isempty(temp))
                [sdf{e, c}(t,:), t_] = spikeTrainDensity(-considerWindowMs/1000, considerWindowMs/1000, temp/(sampleRate), 0.3, 5);
%             figure(c);
%             hold on; plot(temp, t, '*b');
            else
                sdf{e,c}(t,:) = zeros(1, 3532); %spikeTrainDensity(-considerWindowMs/1000, considerWindowMs/1000, 0.0, 1, 20);
            end
        end
%         axis([-considerWindowMs considerWindowMs 0 20])
    end
end

%% Plot raster and SDF for each electrode/recording relation

close all

mkdir([pathDir '/mat_figs']);
mkdir([pathDir '/jpgs']);

chPerPlot = 256;

factor = 5000;

for e = 1:256 %ChannelsToConsider
%     image = squeeze(binFrames(e, :,:))';
%     image = image(end:-1:1,:,:);
%     image = image(:,end:-1:1,:);
e
    n=0;
    chans2cons = [];
    for ch = ChannelsToConsider
        toPlot = (sdf{e,ch});
        toPlot = mean(toPlot);
        if(~isempty(toPlot) && ( max(toPlot(floor(end/2:end))) - max(toPlot(floor(1:end/2))) ) > 5e-2)
            n=n+1;
            chans2cons = [chans2cons ch];            
        end
    end


    if(n>0)

        c=0;
        gcf = figure(e);  
        set(gcf, 'Position', [0 0 1000 1000]);
        set(gcf, 'Visible', 'off');

        for ch = chans2cons
            c=c+1;
            toPlot = (sdf{e,ch});

	%%% Plots raster and SDF
                subplot(n, 2, 2*c);%-chPerPlot*floor(c/chPerPlot) 

                    for t=1:maxTrials(e)
                        if(~isempty(rasterData{e, ch}{t}))
                            hold on; plot(rasterData{e, ch}{t}/sampleRate, t, '.b');
                            axis([-considerWindowMs/1000 considerWindowMs/1000 0 maxTrials(e)])
                        end 
                    end
                    title(num2str(ch));

                    subplot(n, 2, 2*c-1); plot(t_, mean(sdf{e,ch}));
                    title(num2str(ch));
                                    
    %                 axis([-considerWindowMs considerWindowMs 0 20])
        end
            suptitle(['Electrode stimulated:' num2str(e)]);
            drawnow;
            
            saveas(gcf, [pathDir 'jpgs/SpikeDensityFunction_channel_' num2str(e) '_raster'], 'jpeg')
            saveas(gcf, [pathDir 'mat_figs/SpikeDensityFunction_channel_' num2str(e) '_raster'], 'fig')

        close all;
    end
end
% close all;


%% Plot images with electrode lication and corresponding activity
close all

plotDelays = 1;
    
mkdir([pathDir '/mat_figs']);
mkdir([pathDir '/jpgs']);

chPerPlot = 256;

factor = 5000;

image = imread([pathMain '2016_04_08_BB161i_OG_tdt_4x_4sec_13.6x.tif']);
image = squeeze(image(:,:,1))'/1;

for e = 1:256 %ChannelsToConsider

    n=0;
    chans2cons = [];
    for ch = ChannelsToConsider
        toPlot = (sdf{e,ch});
        toPlot = mean(toPlot);
        if(~isempty(toPlot) && ( max(toPlot(floor(end/2:end))) - max(toPlot(floor(1:end/2))) ) > 5e-2)
            n=n+1;
            chans2cons = [chans2cons ch];            
        end
    end
[e n]

    if(n>0)

        %%% Plot images to show chan/elec positions
        
        gcf = figure(e); 
        set(gcf, 'Visible', 'off');
        imshow(image);  

        hold on; scatter(electrode(find(MEA_MAP1 == electrode(e).ID)).coord_TIFF(1)-15, electrode(find(MEA_MAP1 == electrode(e).ID)).coord_TIFF(2)+5, 80, 'y', 'Filled');

%             for chan = 1:256
%                 k = find(MEA_MAP1 == electrode(chan).ID);
%                 hold on; scatter(electrode(k).coord_TIFF(1), electrode(k).coord_TIFF(2), 5, 'r');
%             end
%             drawnow;
            

    %     figure(e);
    c=0;
        for ch = chans2cons
    %             figure( ceil(c/chPerPlot) );

            c=c+1;
            toPlot = (sdf{e,ch});

            impData = mean(sdf{e,ch});
            impData = impData(floor(end/2:end));

            [mx, ~] = max(impData);
            [indmx] = min(find(impData == max(impData)));
            t_(indmx+floor(end/2))

            k = find(MEA_MAP1 == electrode(ch).ID);                
            if(plotDelays)
                hold on; hs = scatter(electrode(k).coord_TIFF(1)-15, electrode(k).coord_TIFF(2)+5, t_(indmx+floor(end/2))*20*factor, 'r', 'Filled');
                text(electrode(k).coord_TIFF(1)-15, electrode(k).coord_TIFF(2)+5, num2str(1000*t_(indmx+floor(end/2))), 'Color', cols(10,:), 'FontWeight', 'Bold', 'FontAngle', 'italic');
            else
                hold on; hs = scatter(electrode(k).coord_TIFF(1)-15, electrode(k).coord_TIFF(2)+5, mx*factor, 'b', 'Filled');
            end
            alpha(hs, 0.2);
        end
        suptitle(['Electrode stimulated:' num2str(e)]);

        if(plotDelays)
            saveas(gcf, [pathDir '/jpgs/SpikeDensityFunction_channel_image_flors_delay_' num2str(e) ], 'jpeg')
            saveas(gcf, [pathDir '/mat_figs/SpikeDensityFunction_channel_image_flors_delay_' num2str(e) ], 'fig')
        else
            saveas(gcf, [pathDir '/jpgs/SpikeDensityFunction_channel_image_flors_' num2str(e) ], 'jpeg')
            saveas(gcf, [pathDir '/mat_figs/SpikeDensityFunction_channel_image_flors_' num2str(e) ], 'fig')
        end
        close all;
    end
end
% close all;


%% Plot all elecs in one image

close all

plotDelays = 1;

mkdir([pathDir '/mat_figs']);
mkdir([pathDir '/jpgs']);

chPerPlot = 256;

factor = 5000;

[image, ~] = imread([pathMain '2016_04_08_BB161i_OG_tdt_4x_4sec_13.6x.tif']);
% image = squeeze(image(:,:,1))'/1;
image1(:,:, 1) = squeeze(image(:,:,1))';
% image1(:,:, 2) = squeeze(image(:,:,2))';
% image1(:,:, 3) = squeeze(image(:,:,3))';

gcf = figure(1); 
imshow(image1);

cols = hsv(38);

colCount = 0;            
theta_rot = 0*pi/180;
chansLink = zeros(1,256);

for e = 1:256 %ChannelsToConsider

    e
    n=0;
    chans2cons = [];
    for ch = ChannelsToConsider
        toPlot = (sdf{e,ch});
        toPlot = mean(toPlot);
        if(~isempty(toPlot) && ( max(toPlot(floor(end/2:end))) - max(toPlot(floor(1:end/2))) ) > 5e-2)
            n=n+1;
            chans2cons = [chans2cons ch];
            chansLink(ch) = chansLink(ch) + 1;
        end
    end
    
    if(n>0)

        colCount = colCount+1;
        %%% Plot images to show chan/elec positions
        
        k = find(MEA_MAP1 == electrode(e).ID);
        pointElectrode = [electrode(k).coord_TIFF(1)*cos(theta_rot) + electrode(k).coord_TIFF(2)*sin(-theta_rot)-25, electrode(k).coord_TIFF(1)*sin(theta_rot) + electrode(k).coord_TIFF(2)*cos(theta_rot)+5];
        
        hold on; scatter(pointElectrode(1), pointElectrode(2), 50, 'MarkerFaceColor', cols(colCount,:), 'MarkerEdgeColor', cols(colCount,:));
        c=0;
        for ch = chans2cons
            c=c+1;
            toPlot = (sdf{e,ch});

            impData = mean(toPlot);
            impData = impData(floor(end/2:end));
            [mx, ~] = max(impData);
            [indmx] = min(find(impData >= max(impData)/2));
            t_(indmx+floor(end/2))
            
            k = find(MEA_MAP1 == electrode(ch).ID);           
            pointChannel = [electrode(k).coord_TIFF(1)*cos(theta_rot) + electrode(k).coord_TIFF(2)*sin(-theta_rot)-25, electrode(k).coord_TIFF(1)*sin(theta_rot) + electrode(k).coord_TIFF(2)*cos(theta_rot)+5];
            hold on; 
%             hs = scatter(pointChannel(1), pointChannel(2), max(impData)*factor, 'MarkerEdgeColor', cols(colCount,:), 'LineWidth', 1.5);
            if(plotDelays)
                hs = scatter(pointChannel(1), pointChannel(2), (t_(indmx+floor(end/2)))*10*factor, 'MarkerEdgeColor', cols(colCount,:));
%                 line([pointElectrode(1) pointChannel(1)], [pointElectrode(2) pointChannel(2)], 'Color', cols(colCount,:), 'LineWidth',  t_(indmx+floor(end/2))*100);
%                 text(pointChannel(1)+(chansLink(ch)-1)*50, pointChannel(2), num2str(1000*t_(indmx+floor(end/2))), 'Color', cols(colCount,:), 'FontWeight', 'Bold', 'FontAngle', 'italic');
            else
                hs = scatter(pointChannel(1), pointChannel(2), (mx)*factor, 'Filled', 'MarkerEdgeColor', cols(colCount,:), 'MarkerFaceColor', cols(colCount,:));
            end
            alpha(hs, 0.2);
            
        end
        drawnow;
%         suptitle(['Electrode stimulated:' num2str(e)]);

%         close all;
    end
end

if(plotDelays)
    saveas(gcf, [pathDir '/jpgs/SpikeDensityFunction_channel_image_flors_all_delays'], 'jpeg')
    saveas(gcf, [pathDir '/mat_figs/SpikeDensityFunction_channel_image_flors_all_delays'], 'fig')
else
    saveas(gcf, [pathDir '/jpgs/SpikeDensityFunction_channel_image_flors_all'], 'jpeg')
    saveas(gcf, [pathDir '/mat_figs/SpikeDensityFunction_channel_image_flors_all'], 'fig')
end

close all;

%% Save workspace

save(['C:\Users\himanshu-idv\POST_DOC\DATA\MEA_DATA\Primates\08_04_2016\SpikeTimes\2016_04_20_BA901h_OD1_gsmmix_DMD600i100_10ms' num2str(spotSize) 'um_' num2str(timeStim) 'ms_' num2str(freq) 'hz\analyseElectrodeWise.mat']);

exist(['C:\Users\himanshu-idv\POST_DOC\DATA\MEA_DATA\Primates\08_04_2016\SpikeTimes\2016_04_20_BA901h_OD1_gsmmix_DMD600i100_10ms' num2str(spotSize) 'um_' num2str(timeStim) 'ms_' num2str(freq) 'hz\analyseElectrodeWise.mat'])
%             line([pointElectrode(1) pointChannel(1)], [pointElectrode(2) pointChannel(2)], 'Color', cols(colCount,:), 'LineWidth',  max(impData)*5);
