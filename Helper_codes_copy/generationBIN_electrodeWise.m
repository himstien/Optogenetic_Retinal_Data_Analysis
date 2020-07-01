function [nomFichier, binFrames] = generationBIN_electrodeWise(electrode, distanceBetweenElectrodes, DMD, fileName, c)

%%V1 OK?? . Generate electrode wise stimulation given DMD coordinates

% removed video making and moved to vecs

%% intialisation du BIN
%Tailles du DMD

factor = 1;

screenWidth=1024;
screenHeight=768;

numElectrodes = size(electrode,2);

numFrames = numElectrodes;

nomFichier=[fileName '_MEA_' num2str(2*2.27*distanceBetweenElectrodes/factor) 'um_' int2str(c(3)) '.' int2str(c(2)) '.' int2str(c(1)) '_' int2str(c(4)) 'h' int2str(c(5)) '.bin'];
fid=fopen(nomFichier,'w','ieee-le'); %Ouverture du nouveau fichier BIN
%fid=fopen('C:\Program Files\ALP-4.1\ALP-4.1 high-speed% API\Samples\BINVECS\Kevin\Bin\STA_pixels_eclairesTEST2.bin','w','ieee-le');

disp(['Writing BIN file to: ' nomFichier]);

StimData=int16([screenWidth screenHeight numFrames 1]); %[screensizeX screensizeY numFrames 1]
fwrite(fid,StimData,'int16');

bitf = bitframe(fid);



Nx = 768; Ny = 1024;
[x, y] = meshgrid(1:Nx, 1:Ny);


%% Recup�ration de la matrice souhait�e

Frame=uint8(0*ones(screenWidth,screenHeight));
bitf.add(uint8(Frame)); % Ecriture de la frame pour bitcompression


%% Create frames for the bin

h = waitbar(0,'Writing files...');

%% Debug
%     Frame=(zeros(screenWidth,screenHeight));
% 
% for i=1:numFrames
%     
%     MEA2Stim = 255*(hypot(x-(electrode(i).coord_DMD(2)), y-(electrode(i).coord_DMD(1)))< distanceBetweenElectrodes/5);
%     Frame= Frame + MEA2Stim;    %
% end
% 
% % DMD2STIM = 255*(hypot(x-DMD.coord.Left(2), y-DMD.coord.Left(1))< distanceBetweenElectrodes/2);
% %     Frame = Frame + DMD2STIM;
% %     
% %     DMD2STIM = 255*(hypot(x-DMD.coord.Right(2), y-DMD.coord.Right(1))< distanceBetweenElectrodes/2);
% %     Frame = Frame + DMD2STIM;
% %     
% %     DMD2STIM = 255*(hypot(x-DMD.coord.Bottom(2), y-DMD.coord.Bottom(1))< distanceBetweenElectrodes/2);
% %     Frame = Frame + DMD2STIM;    
% 
%     Frame = Frame(end:-1:1,:); %Sym�trie horyzontale des frames
%     Frame = Frame(:, end:-1:1); %Sym�trie horyzontale des frames
% 
%     figure; imshow(Frame);

%% 
for frame = 1:numFrames
    tic
    Frame=(zeros(screenWidth,screenHeight));

%     DMD2STIM = 255*(hypot(x-DMD.coord.Left(2), y-DMD.coord.Left(1))< distanceBetweenElectrodes/2);
    
    MEA2Stim = 255*(hypot(x-(electrode(frame).coord_DMD(2)), y-(electrode(frame).coord_DMD(1)))< distanceBetweenElectrodes/factor);
    Frame= Frame + MEA2Stim;    %
    
    Frame = Frame(end:-1:1,:); %Sym�trie horyzontale des frames
    Frame = Frame(:, end:-1:1); %Sym�trie verticale des frames
% % % % % % % % % %      Frame=abs((255*ones(screenWidth,screenHeight))-Frame);
    bitf.add(uint8(Frame)); %�criture de la frame pour bitcompression
    
    t = toc;
%     disp(['Time remaining: ' num2str( (numFrames - frame)*t) ' sec.']);
    perc = (frame*100)/numFrames;
    waitbar(perc/100, h, sprintf('Writing bin file %d%%', uint8(perc)));
    binFrames(frame, :, :) = Frame;
end

fclose(fid); % fermeture du fichier BIN
close(h);

