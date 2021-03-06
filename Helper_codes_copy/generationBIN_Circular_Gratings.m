% function  = generationBIN_Gratings()


%% intialisation du BIN
%Tailles du DMD
clear;

close all;

converstionFactor = 2.27; % 1 pixel = 2.27 um

DMD = 1;

if(DMD==1)
    screenWidth=1024;
    screenHeight=768;
else
    screenWidth = 1920;
    screenHeight = 1080;
end

[x, y] = meshgrid(1:screenHeight, 1:screenWidth);

frameReference = [0 0 0 0];

% gcf = figure;

count = 0;

date = '07_10_2016';

disp('Creating frames..');

c=0;
for f = [4 6 8 12 16 32]
    gratingSize = round(screenHeight/(2*f));
    c = c+gratingSize;
end

Frames = zeros(c*2, screenWidth, screenHeight);

useFoveaCenter = 1;
if (useFoveaCenter)
    defaultPath = ['/home/icub/POST_DOC/DATA/' date '_opto_primate'];
    filePath = [defaultPath '/'];
    test = input('Enter test number for calibration figures: ', 's');
    fileName = 'test';

    DMD_fileName = [filePath fileName test 'c.tif'];
    [X_DMD, Y_DMD, DMD_Image] = calibrationRepere_HA(DMD_fileName);

    floresence_fileName = [filePath fileName test 'f.tif'];

    center = getFoveaCenter(floresence_fileName, X_DMD, Y_DMD);
else
    center = [screenHeight/2 screenWidth/2];
    center = [300 668];
end

center
%%
for factor = [4 6 8 12 16 32] %[2 4 8 10 32]
    
    factor
    gratingSize = round(screenHeight/(2*factor));
    gratingSize_um = gratingSize*(converstionFactor); %um

    % make 1st frame

    for phase = 1:gratingSize
        range = gratingSize+phase-1:gratingSize*2:screenHeight-1;

        frame = zeros(screenWidth, screenHeight);
        
        for r = 1:numel(range)
            circle = 255* (hypot(x-center(2), y-center(1)) > range(r)-gratingSize & hypot(x-center(2), y-center(1)) < range(r));
%             inCircle = 255* (hypot(x-center(1), y-center(2)) < range(r));

            frame = frame + circle;       
        end        


        count = count+1;

        Frames(count, :, :) = frame;
        frameReference = [frameReference; [gratingSize gratingSize_um phase 1] ];

%         imshow(frame);
%         drawnow;

        count=count+1;
        Frames(count, :, :) = 255 - frame;
        frameReference = [frameReference; [gratingSize gratingSize_um phase 0] ];

%         imshow(255 - frame);
%         drawnow;

    end
end

numFrames = count;
disp(['Created ' num2str(count) ' frames to write']);

%%
nomFichier=['Flipping_Circular_Gratings_DMD_centered_Fovea_' num2str(screenWidth) '_' num2str(screenHeight) '_' date '.bin'];
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

    Frame = Frame(end:-1:1,:); %Sym�trie horyzontale des frames
    Frame = Frame(:, end:-1:1); %Sym�trie verticale des frames
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

save(['Circular_Grating_Stimulation_binFile_Reference_normal_' date '.mat'], 'frameReference', '-v7.3');
% save('Circular_Grating_Stimulation_binFile_Frames_normal.mat', 'Frames', '-v7.3');

%%

I = imread(floresence_fileName);
imshow(squeeze(I(:,:,1))');
hold on; scatter(center(1), center(2), 30, 'o', 'Filled', 'MarkerFaceColor', [1 0 0]);


%% CREATE VEC

periode = 1;
nbReps = 50;

% Total de frames diffus�s
% nFrames = size(binFrames,1);

totalnFrames = nbReps*(numFrames*periode);

frs = 1:numFrames;

% nomFichier=['oriented_movingBars_reps_' num2str(nbReps) '_timeOn_' num2str(periode) '.vec'];

nomFichier=['circular_Grating_' num2str(nbReps) '_timeOn_' num2str(periode) '_' date '.vec'];


fid=fopen(nomFichier,'w');
OutputVec=[0 totalnFrames 0 0 0];
fprintf(fid,'%g %g %g %g %g\n',OutputVec);

disp(['Writing vec file to location: ' nomFichier]);

Nx = 1024; Ny = 768;

[x, y] = meshgrid(1:Nx, 1:Ny);
% Frames=zeros(screenWidth,screenHeight,numFrames+1); % on cr�� la matrice 3D qui va contenir toutes les frames

%%%% MAIN LOOP
%baseline de 'periode' ms


for reps = 1:nbReps

    indices = 1:numFrames;
    for i=1:10
        OutputVec=[0 0 0 0 0];
        fprintf(fid,'%g %g %g %g %g\n',OutputVec);
    end

    for frame=1:numFrames        
        for stim=1:periode
            a = frs(indices(frame));
            OutputVec=[0 a 0 0 0];
            fprintf(fid,'%g %g %g %g %g\n',OutputVec);
        end
    end
end

fclose(fid);

disp('Finished writing vec file...');
    