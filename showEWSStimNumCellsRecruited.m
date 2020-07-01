% close all
% for e = 1:256
% for c = ChannelsToConsider
% if( (max(squeeze(sdf(e,c,100:500)) ) - mean(squeeze(sdf(e,c,100:500)) ))) > 300
% figure(e);
% hold on;
% plot(squeeze(sdf(e,c, 100:500)));
% end
% end
% 
% end


%%


size_ = {'50um', '75um', '100um'};

tempDir = '/media/icub/B4A68AB8A68A7B1E/icub/POST_DOC/DATA/20_04_2016_opto_primate/SpikeTimes/';

image = imread(['./' 'test6f.tif']);
image = squeeze(image(:,:,1))'/1;

colmaps = jet(256);
numRecruits = [];

% gcf = figure(size__); 
% set(gcf, 'Visible', 'on');
% imshow(image');
%             
% for e = [4:252] %ChannelsToConsider   
%         hold on;
%         text(electrode(find(MEA_MAP1 == electrode(e).ID)).coord_TIFF(2)-55, electrode(find(MEA_MAP1 == electrode(e).ID)).coord_TIFF(1)-55, num2str(e ), 'Color', 'y');
% 
% end
%%

close all;

for size__ = 3

    tempDir_ = ['2016_04_20_BA901h_OD1_gsmmix_DMD600i100_EWS_' size_{size__} '_25Hz/'];

    load(['SpikeTimes/' tempDir_ '/analyseElectrodeWise.mat']);
   %Plot images with electrode lication and corresponding activity
%     close all

    plotDelays = 1;

    % mkdir([pathDir '/mat_figs']);
    % mkdir([pathDir '/jpgs']);

    chPerPlot = 256;

    factor = 5000;

    image = imread(['./' 'test6f.tif']);
    image = squeeze(image(:,:,1))'/1;

    colmaps = jet(256);
    numRecruits = [];

    gcf = figure(size__); 
            set(gcf, 'Visible', 'on');
            imshow(image');  

drawnow();

    for e = [59:61 73:79 88:95 103:110 119:125 135:140 151:156 167:172 183:189 200:205 216:221 234:237 251 252] %ChannelsToConsider

%         if( isempty( find(e == [59 60 61 216 251 252 250 17 20 233 55 71 82 193 190 192 46 173 214 231 167  49 51 44 130 68 100 240 158 103 227 14 213 248 119 135 151 16 24 27 31 47 144 176 165 118 166 142 64 63 95 57 73 33 162  117 242 177]) ))
            n=0;
            chans2cons = [];
            numRecruits(e) = 0;
            for ch = [59:61 73:79 88:95 103:110 119:125 135:140 151:156 167:172 183:189 200:205 216:221 234:237 251 252]%ChannelsToConsider
                
                toPlot = squeeze(sdf(e,ch,:));
        %         toPlot = mean(toPlot);
                if(~isempty(toPlot) && ( max(toPlot(floor(end/2:end))) - max(toPlot(floor(1:end/2))) ) > 250)
                    n=n+1;
                    chans2cons = [chans2cons ch];
                    numRecruits(e) = n;
                    
%                     figure(10*size__);
%                     hold on; plot(toPlot);
                end

            end
        % [e n]

            if(numRecruits(e)>0)

                %%% Plot images to show chan/elec positions


        %         hold on; scatter(electrode(find(MEA_MAP1 == electrode(e).ID)).coord_TIFF(1)-15, electrode(find(MEA_MAP1 == electrode(e).ID)).coord_TIFF(2)+5, 80, 'y', 'Filled');


                hold on; hs = scatter(electrode(find(MEA_MAP1 == electrode(e).ID)).coord_TIFF(2)-45, electrode(find(MEA_MAP1 == electrode(e).ID)).coord_TIFF(1)-55, 100*(numRecruits(e)), colmaps(e,:), 'Filled', 'MarkerFaceAlpha', 1 );
                text(electrode(find(MEA_MAP1 == electrode(e).ID)).coord_TIFF(2)-50, electrode(find(MEA_MAP1 == electrode(e).ID)).coord_TIFF(1)-55, num2str(numRecruits(e)), 'Color', 'k' );
%                 for nnn = 1:numel(chans2cons)
%                     scatter(electrode(find(MEA_MAP1 == electrode(chans2cons(nnn)).ID)).coord_TIFF(2)-45, electrode(find(MEA_MAP1 == electrode(chans2cons(nnn)).ID)).coord_TIFF(1)-55, 150*numel(numRecruits(e)), colmaps(e,:), 'MarkerFaceAlpha', 1 );
%                 end

%                 text(electrode(find(MEA_MAP1 == electrode(e).ID)).coord_TIFF(2)-55, electrode(find(MEA_MAP1 == electrode(e).ID)).coord_TIFF(1)-55, num2str(numRecruits(e) ), 'Color', 'k');
%                 
%                 text(electrode(find(MEA_MAP1 == electrode(e).ID)).coord_TIFF(2)-55, electrode(find(MEA_MAP1 == electrode(e).ID)).coord_TIFF(1)-55, num2str(e ), 'Color', 'm');
drawnow;
        %             for chan = 1:256
        %                 k = find(MEA_MAP1 == electrode(chan).ID);
        %                 hold on; scatter(electrode(k).coord_TIFF(1), electrode(k).coord_TIFF(2), 5, 'r');
        %             end
        %             drawnow;


            %     figure(e);
        %     c=0;
        %         for ch = chans2cons
        %     %             figure( ceil(c/chPerPlot) );
        % 
        %             c=c+1;
        %             toPlot = (sdf{e,ch});
        % 
        %             impData = mean(sdf{e,ch});
        %             impData = impData(floor(end/2:end));
        % 
        %             [mx, ~] = max(impData);
        %             [indmx] = min(find(impData == max(impData)));
        %             t_(indmx+floor(end/2))
        % 
        %             k = find(MEA_MAP1 == electrode(ch).ID);                
        %             if(plotDelays)
        %                 hold on; hs = scatter(electrode(k).coord_TIFF(1)-15, electrode(k).coord_TIFF(2)+5, t_(indmx+floor(end/2))*20*factor, 'r', 'Filled');
        %                 text(electrode(k).coord_TIFF(1)-15, electrode(k).coord_TIFF(2)+5, num2str(1000*t_(indmx+floor(end/2))), 'Color', cols(10,:), 'FontWeight', 'Bold', 'FontAngle', 'italic');
        %             else
        %                 hold on; hs = scatter(electrode(k).coord_TIFF(1)-15, electrode(k).coord_TIFF(2)+5, mx*factor, 'b', 'Filled');
        %             end
        %             alpha(hs, 0.2);
        %         end
        %         suptitle(['Electrode stimulated:' num2str(e)]);
        % 
        %         if(plotDelays)
        %             saveas(gcf, [pathDir '/jpgs/SpikeDensityFunction_channel_image_flors_delay_' num2str(e) ], 'jpeg')
        %             saveas(gcf, [pathDir '/mat_figs/SpikeDensityFunction_channel_image_flors_delay_' num2str(e) ], 'fig')
        %         else
        %             saveas(gcf, [pathDir '/jpgs/SpikeDensityFunction_channel_image_flors_' num2str(e) ], 'jpeg')
        %             saveas(gcf, [pathDir '/mat_figs/SpikeDensityFunction_channel_image_flors_' num2str(e) ], 'fig')
        %         end
        %         close all;
            end
        end
    
    % close all;
end
