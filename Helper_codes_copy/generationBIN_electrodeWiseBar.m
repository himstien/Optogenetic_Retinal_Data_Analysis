function   generationBIN_electrodeWiseBar(widthBar, angles, direction)

%%V1 OK?? . Generate electrode wise stimulation given DMD coordinates

% removed video making and moved to vecs

%% intialisation du BIN
%Tailles du DMD

close all;
% clear;

writeFile = 0;

% widthBar = 100;
factor = 1;

screenWidth=1024;
screenHeight=768;

numFrames = screenHeight;

reverse = strcmp(direction, 'reverse');

if(writeFile)
    
    if(reverse)
        nomFichier=['Oriented_MOVING_BAR_' num2str(widthBar*2.27) '_um_reverse' '.bin'];
    else
        nomFichier=['Oriented_MOVING_BAR_' num2str(widthBar*2.27) '_um_' '.bin'];
    end
    
    fid=fopen(nomFichier,'w','ieee-le'); %Ouverture du nouveau fichier BIN
    %fid=fopen('C:\Program Files\ALP-4.1\ALP-4.1 high-speed% API\Samples\BINVECS\Kevin\Bin\STA_pixels_eclairesTEST2.bin','w','ieee-le');

    disp(['Writing BIN file to: ' nomFichier]);

    StimData=int16([screenWidth screenHeight numFrames 1]); %[screensizeX screensizeY numFrames 1]
    fwrite(fid,StimData,'int16');

    bitf = bitframe(fid);
    
    Frame=uint8(0*ones(screenWidth,screenHeight));
    bitf.add(uint8(Frame)); % Ecriture de la frame pour bitcompression

end

Nx = 768; Ny = 1024;
[x, y] = meshgrid(1:Nx, 1:Ny);



%% Create frames for the bin

h = waitbar(0,'Writing files...');

center = [screenWidth/2 screenHeight/2];

%% 

% binFrames = zeros(numel(angles)*numFrames, screenWidth, screenHeight);
frameReference = [0 0 0 0];

% jump = 4;
c=0;
for ang = angles
    for frame = 1:numFrames
        tic
        c=c+1;
        
        Frame=(zeros(screenWidth,screenHeight));

    %     DMD2STIM = 255*(hypot(x-DMD.coord.Left(2), y-DMD.coord.Left(1))< distanceBetweenElectrodes/2);
        toWrite = round(frame:frame+widthBar/factor);
        toWrite = toWrite(toWrite>0 & toWrite <= screenHeight);
        Frame(:,toWrite) = 255;
        Frame = imrotate(Frame, ang, 'crop');

    %     Frame = Frame( size(Frame,1)-screenWidth/2:size(Frame,1)+screenWidth/2-1, size(Frame,2)-screenHeight/2:size(Frame,2)+screenHeight/2-1);

    %     close all
    %     MEA2Stim = 255*(hypot(x-(electrode(frame).coord_DMD(2)), y-(electrode(frame).coord_DMD(1)))< distanceBetweenElectrodes/factor);
    %     Frame= Frame + MEA2Stim;    %

    if(reverse)
        Frame = Frame(end:-1:1,:); %Sym�trie horyzontale des frames
        Frame = Frame(:, end:-1:1); %Sym�trie verticale des frames
    end        

%     imshow(Frame);
%     drawnow;
%     pause(0.1);
    
    % % % % % % % % % %      Frame=abs((255*ones(screenWidth,screenHeight))-Frame);

        if(writeFile)
            bitf.add(uint8(Frame)); %�criture de la frame pour bitcompression
        end
        
        frameReference = [frameReference; [widthBar widthBar*2.27 ang frame] ];

        t = toc;
    %     disp(['Time remaining: ' num2str( (numFrames - frame)*t) ' sec.']);
        perc = (c*100)/(numel(angles)*numFrames);
        waitbar((perc/100), h, sprintf('Writing bin file %d%%', uint8(perc)));
%         binFrames(c, :, :) = Frame;
    end
end

    if(writeFile)
        fclose(fid); % fermeture du fichier BIN
    end

    
close(h);

%%
% save([nomFichier(1:end-4) '_frames.mat'], 'binFrames', '-v7.3');
if(reverse)
    save(['Oriented Bar_binFile_Reference_reverse.mat']);
end
