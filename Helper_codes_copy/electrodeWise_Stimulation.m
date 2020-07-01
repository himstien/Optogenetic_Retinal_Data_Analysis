%Create 3-D grid of implant pixels from p.tiff image
%% Initialize
close all;
clear variables;
clear global;
clc;

%date= 'test'; %input('Mettre la date de l''exp�rience (ex 20150120): ','s');

% MAXPIXELSON = min(MAXPIXELSON, MAXPIXELSON_Vec); % use 

defaultPath = 'C:\Users\EQS8\Desktop\Test_greg';
filePath = [defaultPath '\'];
test = input('Enter test number for calibration figures: ', 's');
fileName = 'test';

createVec = true;
createBin = true;
createVideo = false;
saveWS = true;

%% Get file name
% [fileName, filePath, s] = uigetfile('*.tiff');
% indexTemp = findstr(fileName, '.tiff');
% fileName = fileName(1:indexTemp-2);
% clear indexTemp
%% Get coordinates of DMD from image

DMD_fileName = [filePath fileName test 'c.tif'];
[X_DMD, Y_DMD, DMD_Image] = calibrationRepere_HA(DMD_fileName);
DMD.coord.topLeft = [X_DMD(1) Y_DMD(1)];
DMD.coord.bottomLeft = [X_DMD(2) Y_DMD(2)];
DMD.coord.topRight = [X_DMD(3) Y_DMD(3)];

DMD.coord.Left = coordTIFF2DMD_v2(X_DMD, Y_DMD, DMD.coord.topLeft', 2);
DMD.coord.Bottom = coordTIFF2DMD_v2(X_DMD, Y_DMD, DMD.coord.bottomLeft', 2);
DMD.coord.Right= coordTIFF2DMD_v2(X_DMD, Y_DMD, DMD.coord.topRight', 2);

%% Get electrode locations

pictureNameCalibrationElectrods =  [filePath fileName test 'e.tif'];
[electrode, distanceBetweenElectrodes, electrode_Image] = electrodes(X_DMD, Y_DMD, pictureNameCalibrationElectrods);

%% Generate frames and bin file 
% And generates VEC file for Gregory


stim_fileName = './'; %[filePath fileName];
bin_fileName = [stim_fileName fileName test];
vec_fileName = [stim_fileName fileName test];

c=clock;

if(createBin)
    [~, binFrames] = generationBIN_electrodeWise(electrode, distanceBetweenElectrodes, DMD, bin_fileName,c);
end


if(createVec)
    tps_stim=0;
    while (tps_stim < 2 || mod(tps_stim,2)~=0)
        tps_stim=round(input('Temps de stimulation en ms (min=2ms, nombre pair) : '));
    end
%     tps_stim=tps_stim; %Car diffusion � 500Hz
    
    frequenceStim=0;
    while (frequenceStim<1)
        frequenceStim=round(input('Fr�quence de stimulation en Hz : '));
    end
    periode=1000/frequenceStim;% Car diffusion des frames � 500Hz
    
%     nReps = input('Entrer number de repetition: ');
    
    

%     [nomFichierVEC] = generationVEC_STA_multiPixel(nReps,periode,tps_stim,implant, 402, c, vec_fileName);    
      nReps = input(['Entrer nombre de repetition pour stimulation : ']);
      [nomFichierVEC] = generationVEC_electrodeWise(nReps, periode, tps_stim, electrode, distanceBetweenElectrodes, date, c, vec_fileName, binFrames, createVideo);
end


%% PLOTs

hold off;
% imshow(DMD_Image/4 );
imshow(DMD_Image/4 + electrode_Image/2);
hold on;

ang=0:0.01:2*pi; 

drawCircle([DMD.coord.topLeft(1) DMD.coord.topLeft(2)], 3, 'r');
drawCircle([DMD.coord.bottomLeft(1), DMD.coord.bottomLeft(2)], 3, 'g');
drawCircle([DMD.coord.topRight(1), DMD.coord.topRight(2)], 3, 'b');

% drawCircle([DMD.coord.Left(1), DMD.coord.Left(2)], 3, 'm');
% drawCircle([DMD.coord.Bottom(1), DMD.coord.Bottom(2)], 3, 'c');
% drawCircle([DMD.coord.Right(1), DMD.coord.Right(2)], 3, 'w');

for k=1:256
    drawCircle([electrode(k).coord_TIFF(1),electrode(k).coord_TIFF(2)], 3, 'k');
%     drawCircle([electrode(k).coord_DMD(1),electrode(k).coord_DMD(2)], 3, 'y');
end
% for n = 1:size(implant.matrixStimPairsImage,1)
%     figure(2);
%     subplot(ceil(sqrt(size(implant.matrixStimPairsImage,1))), ceil(sqrt(size(implant.matrixStimPairsImage,1))), n);
%     imshow(pixelsImplant);
%     hold on;    
%     pause(0.20)
% %     drawCircle([implant.matrixStimPairsImage{1}(1,1), implant.matrixStimPairsImage{1}(1,2)], implant.distBetweenPixels*2/2 , 'y');
%     for s = 1:size(implant.matrixStimPairsImage{n},1)                
%         drawCircle([implant.matrixStimPairsImage{n}(s,1), implant.matrixStimPairsImage{n}(s,2)], implant.distBetweenPixels/2 , 'r');
%     end
% end
saveas(gcf, ['_final.tiff']);
saveas(gcf, ['_final.fig']);

%% 
if(saveWS)
    nomWorkSpace=['./', fileName, test, '_',date,'_', num2str(tps_stim*2),'ms_' , int2str(c(4)),'h',int2str(c(5)) ];
    save(nomWorkSpace);
    disp(['Sauvegarde des variables OK : ',nomWorkSpace]);
end
%% Temp stuff
