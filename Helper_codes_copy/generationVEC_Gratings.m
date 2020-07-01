
converstionFactor = 2.27; % 1 pixel = 2.27 um

screenWidth=1024;
screenHeight=768;

load('Circular_Grating_Stimulation_binFile_Reference_normal.mat');
frameReference = frameReference(2:end,:);

gratingSize_ = input('Enter grating size in um (nearest set will be selected, enter negative for all):');

if(gratingSize_ > 0)
        allSizes = 0;
        setOfFrames = [];
        for n=1:numel(gratingSize_)
            distances = (frameReference(:,2) - gratingSize_(n)).^2;
            [nearestSize, indexNearest] = min(distances);
            closestPixelSize(n) = frameReference(indexNearest, 1);
            setOfFrames = [setOfFrames; find(frameReference(:,1) == closestPixelSize(n))];
            disp(['Size of grating displayed is: ' num2str(closestPixelSize*converstionFactor) ]);
        end    
else
    allSizes = 1;
    closestPixelSize = frameReference(1,1);
    setOfFrames = 1:size(frameReference, 1);
end

phaseRange = frameReference(setOfFrames, 3);

phaseToShow = input(['Phase to show. Enter a value between ' num2str(min(phaseRange)) ' and ' num2str(max(phaseRange)) '.Enter a negative value to show all.' ]   );


if(phaseToShow < 0)
    if(~allSizes)
        setOfFrames = find(frameReference(:,1) == closestPixelSize);
    end
else
    tempFinds = [];
 	setOfFrames = [];
 
    for n=1:numel(phaseToShow)
        tempFinds = [tempFinds; find(frameReference(:,3) == phaseToShow(n))];
    end        

    if(~allSizes)
        for n = 1:numel(gratingSize_)
            setOfFrames = [setOfFrames; find(frameReference(:,1) == closestPixelSize(n))];
            setOfFrames = intersect( setOfFrames,tempFinds);
        end
    else
        setOfFrames = tempFinds;
    end
end

%%

numFlips = 10;
timeDisplay = 10;
timeRest = 200;

numReps = 20;

prepTime = 100; %ms;

toWrite = [];

% write black frames before start of experiment for prepTime ms

cond = [];
for p = 1:prepTime
    toWrite = [toWrite; 0];
    cond = [cond; 0];
end

% for each repitition
for reps = 1:20
    % for each condition selected above...  The odd frames are the normal
    % once while the next even frame is inverted version
    for condition = 1:2:numel(setOfFrames)
        for flip = 1:numFlips    
            % write normal frame for timeDisplay ms
            for on = 1:timeDisplay
                toWrite = [toWrite; setOfFrames(condition)];
            end
            % write inverted frame for timeDisplay ms            
            for on = 1:timeDisplay
                toWrite = [toWrite; setOfFrames(condition)+1];
            end
        end
            % write black frames to incorporate resting time      
        for rest = 1:timeRest
            toWrite = [toWrite; 0];
        end
    end
end

totalnFrames = numel(toWrite);

nomFichier=['Gratings_' 'size_' num2str(gratingSize_) '_phase_' num2str(phaseToShow(end)) '_reps_' num2str(numReps)  '_timeDisplay_' num2str(timeDisplay)  '_restTime_' num2str(timeRest) '.vec'];

fid=fopen(nomFichier,'w');
OutputVec=[0 totalnFrames 0 0 0];
fprintf(fid,'%g %g %g %g %g\n',OutputVec);

disp(['Writing vec file to location: ' nomFichier]);



%% write frame numbers 


for framesWrite = 1:totalnFrames
	a = toWrite(framesWrite);
	OutputVec=[0 a 0 0 0];
    fprintf(fid,'%g %g %g %g %g\n',OutputVec);    
end

fclose(fid);

disp('Finished writing vec file...');

