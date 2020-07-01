%% G�n�re le vec : FREQUENCE A 500 Hz ! (si pas de division par 10 au d�but)
% nFrames : nombre de frames � diffuser au total
% Nouvelle frame toutes les "periode" ms
function [nomFichier] = generationVEC_electrodeWise(nbReps, periode,tps_stim,date,c, fileName, binFrames, makeVideo)

screenWidth=1024;
screenHeight=768;


% Total de frames diffus�s
nFrames = size(binFrames,1);
frs = 1:nFrames;

nFrames

periode
totalnFrames = (nFrames)*periode*nbReps + periode;
nomFichier=[fileName, '_' date '_tps_' num2str(tps_stim) '_reps_' num2str(nbReps)  '_periode_' num2str(periode)  '_' int2str(c(4)) 'h' int2str(c(5)) '.vec'];

fid=fopen(nomFichier,'w');
OutputVec=[0 totalnFrames 0 0 0];
fprintf(fid,'%g %g %g %g %g\n',OutputVec);

disp(['Writing vec file to location: ' nomFichier]);


Nx = 1080; Ny = 1920;
[x, y] = meshgrid(1:Nx, 1:Ny);
% Frames=zeros(screenWidth,screenHeight,numFrames+1); % on cr�� la matrice 3D qui va contenir toutes les frames

%% Create video

% Create video file to write activity frames
% makeVideo = true;

if(makeVideo)
    vid = VideoWriter([fileName '_MEA_' int2str(c(3)) '.' int2str(c(2)) '.' int2str(c(1)) '_' int2str(c(4)) 'h' int2str(c(5)) '.avi']);
    vid.FrameRate = 30; 
    vid.open();
end


%% main loop

%baseline de 'periode' ms
for i=1:periode

    Frame=(zeros(screenWidth,screenHeight)); % blank frame

    Frame = Frame(:,end:-1:1); %Sym�trie horyzontale des frames
    Frame=abs((255*ones(screenWidth,screenHeight))-Frame); % invert for DMD PHP

    if(makeVideo)    
        f = squeeze(Frame);
        writeVideo(vid, f/255);
    end

    OutputVec=[0 0 0 0 0];
    fprintf(fid,'%g %g %g %g %g\n',OutputVec);
end


for reps = 1:nbReps
%      disp(num2str(reps));
    randomFrameList = randperm(nFrames);
    
    
    for frame=1:nFrames        
        for stim=1:tps_stim
            %         a=find(matrice(frame,:)==1);
%             a = randomFrameList(frame);
            a = frs( randomFrameList(frame));
            
%             Frame=(zeros(screenWidth,screenHeight)); % blank frame
            
%             MEA2Stim = 255*(hypot(x-(electrode(a).coord_DMD(2)), y-(electrode(a).coord_DMD(1)))< distanceBetweenElectrodes/2);
%             Frame= Frame + MEA2Stim;    % Frame with location of the MEA stimulation

%             Frame = Frame(:,end:-1:1); %Sym�trie horyzontale des frames
            
            
            if(makeVideo && reps == 1)  
                Frame = binFrames(a, :,:);
                f = squeeze(Frame);
                writeVideo(vid, f/255);
            end

%             Frame=abs((255*ones(screenWidth,screenHeight))-Frame); % invert

            OutputVec=[0 a 0 0 0];
            fprintf(fid,'%g %g %g %g %g\n',OutputVec);
        end
        
        for black=1:periode-tps_stim
                        
            if(makeVideo)
                Frame=(zeros(screenWidth,screenHeight)); % blank frame
            
                Frame = Frame(:,end:-1:1); %Sym�trie horyzontale des frames
                Frame=abs((255*ones(screenWidth,screenHeight))-Frame); % invert for DMD PHP

                f = squeeze(Frame);
                writeVideo(vid, f/255);
            end
            
            OutputVec=[0 0 0 0 0];
            fprintf(fid,'%g %g %g %g %g\n',OutputVec);
        end
    end
end

fclose(fid);

% for video
if(makeVideo)    
    vid.close();
end


disp('Finished writing vec file...');
    
end

