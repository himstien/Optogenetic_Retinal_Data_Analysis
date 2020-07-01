%% G�n�re le vec : FREQUENCE A 500 Hz ! (si pas de division par 10 au d�but)
% nFrames : nombre de frames � diffuser au total
% Nouvelle frame toutes les "periode" ms

function  generationVEC_electrodeWise_RandomBar(nbReps, periode)

% nReps = 25; %n
% periode = 10; %ms

screenWidth=1024;
screenHeight=768;

frs = [769:769+767];

nFrames = 768;

% Total de frames diffus�s
% nFrames = size(binFrames,1);

totalnFrames = nbReps*(nFrames*periode);

nomFichier=['randomBars_reps_' num2str(nbReps) '_timeOn_' num2str(periode) '.vec'];

fid=fopen(nomFichier,'w');
OutputVec=[0 totalnFrames 0 0 0];
fprintf(fid,'%g %7.f %g %g %g\n',OutputVec);

disp(['Writing vec file to location: ' nomFichier]);

Nx = 1024; Ny = 768;

[x, y] = meshgrid(1:Nx, 1:Ny);
% Frames=zeros(screenWidth,screenHeight,numFrames+1); % on cr�� la matrice 3D qui va contenir toutes les frames


%% main loop

%baseline de 'periode' ms

for reps = 1:nbReps

    indices = randperm(nFrames);
%     for i=1:10
%         OutputVec=[0 0 0 0 0];
%         fprintf(fid,'%g %g %g %g %g\n',OutputVec);
%     end

    for frame=1:nFrames        
        for stim=1:periode
            a = frs(indices(frame));
            OutputVec=[0 a 0 0 0];
            fprintf(fid,'%g %g %g %g %g\n',OutputVec);
        end
    end
end

fclose(fid);

disp('Finished writing vec file...');
    
end

