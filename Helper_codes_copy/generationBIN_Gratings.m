% function  = generationBIN_Gratings()


%% intialisation du BIN
%Tailles du DMD
clear;

close all;

converstionFactor = 2.27; % 1 pixel = 2.27 um

screenWidth=1024;
screenHeight=768;

frameReference = [0 0 0 0];

% gcf = figure;

count = 0;

disp('Creating frames..');

c=0;
for f = [1 2 3 4 6 8 12 16 32]
    gratingSize = round(screenHeight/(2*f));
    c = c+gratingSize;
end

Frames = zeros(c*2, screenWidth, screenHeight);


for factor = [1 2 3 4 6 8 12 16 32] %[2 4 8 10 32]
    
    factor
    gratingSize = round(screenHeight/(2*factor));
    gratingSize_um = gratingSize*(converstionFactor); %um

    % make 1st frame

    for phase = 1:gratingSize
 
        basePattern = zeros(1,gratingSize*2);
        basePattern(1:gratingSize) = 1;
                
        linesOn = repmat(basePattern, 1, screenHeight/(2*gratingSize));
        linesOn_temp = zeros(1, screenHeight+phase-1);
        linesOn_temp(1, phase:end) = linesOn;
        linesOn = linesOn_temp(1:screenHeight);

        frame = repmat(linesOn, screenWidth, 1);    
        frame=uint8(255*frame);
    
        count = count+1;
        Frames(count, :, :) = frame;
        frameReference = [frameReference; [gratingSize gratingSize_um phase 1] ];

        count=count+1;
        Frames(count, :, :) = 255 - frame;
        frameReference = [frameReference; [gratingSize gratingSize_um phase 0] ];
        
    end
%     imshow(frame);
end

numFrames = count;
disp(['Created ' num2str(count) ' frames to write']);

%%
nomFichier=['Flipping_Gratings_DMD_' num2str(screenWidth) '_' num2str(screenHeight) '_.bin'];
fid=fopen(nomFichier,'w','ieee-le'); %Ouverture du nouveau fichier BIN
%fid=fopen('C:\Program Files\ALP-4.1\ALP-4.1 high-speed% API\Samples\BINVECS\Kevin\Bin\STA_pixels_eclairesTEST2.bin','w','ieee-le');

disp(['Writing BIN file to: ' nomFichier]);


StimData=int16([screenWidth screenHeight numFrames 1]); %[screensizeX screensizeY numFrames 1]
fwrite(fid,StimData,'int16');

bitf = bitframe(fid);

Frame=uint8(0*ones(screenWidth,screenHeight));
bitf.add(uint8(Frame)); % Ecriture de la frame pour bitcompression


%% Create frames for the bin

h = waitbar(0,'Writing files...');

 
disp('Writing frames');
for frame = 1:numFrames
    tic
    Frame= squeeze(Frames(frame, :, :));

% 
%     imshow(Frame);
%     pause(0.2);
    
% 	clf(gcf);

    bitf.add(uint8(Frame)); %�criture de la frame pour bitcompression
    
    t = toc;
%     disp(['Time remaining: ' num2str( (numFrames - frame)*t) ' sec.']);
    perc = (frame*100)/numFrames;
    waitbar(perc/100, h, sprintf('Writing bin file %d%%', uint8(perc)));
%     binFrames(frame, :, :) = Frame;
end

fclose(fid); % fermeture du fichier BIN
close(h);

%%
save('Grating_Stimulation_binFile_Reference.mat', 'frameReference', '-v7.3');
save('Grating_Stimulation_binFile_Frames.mat', 'Frames', '-v7.3');
