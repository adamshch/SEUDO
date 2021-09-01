%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LOAD ALLEN DATA INTO SEUDO OBJECT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CHOOSE DATASET

% datInf.Dir   = '/jukebox/tank/jlgauthi/SEUDO/ALLEN_INST/';
datInf.Dir   = '/mnt/bucket/adam/SEUDOdata/ALLEN_INST/';
datInf.Sel   = 4;
datInf.Files = {'511510779_three_session_A_503109347.h5',...
             '511510779_three_session_B_503019786.h5',...  
             '511510779_three_session_C_502634578.h5',...
             '511510848_three_session_A_506823562.h5',...
             '511510848_three_session_B_507691834.h5',...
             '511510848_three_session_C_504809131.h5'};
datInf.DFF   = {'allendata_dff_503109347.h5','allendata_dff_503019786.h5',...
             'allendata_dff_502634578.h5','allendata_dff_506823562.h5',...
             'allendata_dff_507691834.h5','allendata_dff_504809131.h5'};
datInf.ROIs  = {'allendata_rois_503109347.h5','allendata_rois_503019786.h5',...
             'allendata_rois_502634578.h5','allendata_rois_506823562.h5',...
             'allendata_rois_507691834.h5','allendata_rois_504809131.h5'};


fn.Ims = datInf.Files{datInf.Sel};
fn.T   = datInf.DFF{datInf.Sel};
fn.P   = datInf.ROIs{datInf.Sel};

if datInf.Sel==1
    datInf.classFile = '511510779_three_session_A_503109347_SEUDOclassification.mat';
elseif datInf.Sel==4
    datInf.classFile = '511510848_three_session_A_506823562_SEUDOclassification.mat';
%     datInf.classFile = [];
else
    datInf.classFile = [];
end
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GET SIZES

H                = h5info([datInf.Dir fn.Ims],'/data');
datInf.nY        = H.Dataspace.Size(1);
datInf.nX        = H.Dataspace.Size(2);
datInf.nFrames   = H.Dataspace.Size(3);                                    % OLD: nFrames = 10000;
H                = h5info([datInf.Dir fn.T],'/dff');
datInf.nProfiles = H.Dataspace.Size(2);

clear H
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GET MEAN IMAGE

dataForMean = h5read([datInf.Dir fn.Ims],'/data',[1 1 1],[datInf.nY ...
                        datInf.nX 100],[1 1 floor(datInf.nFrames / 100)]); % Get the data necessary to calculate the mean iamge
meanIm      = double(mean(dataForMean,3));                                 % Calculate the mean image

clear dataForMean
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GET TIME COURSES & ROIs

T = h5read([datInf.Dir fn.T],'/dff',[1 1],[datInf.nFrames ...
                                                 datInf.nProfiles],[1 1]); % Extract dF/F time courses
P = double(h5read([datInf.Dir fn.P],'/X',[1 1 1],[datInf.nY datInf.nX ...
                                              datInf.nProfiles],[1 1 1])); % Extract spatial profiles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SET UP FUNCTION

getIm = @(f)double(h5read([datInf.Dir fn.Ims],'/data',[1 1 f],[512 512 1],...
                                                         [1 1 1]))-meanIm; % This function reads in one data frame

if isempty(datInf.classFile)
    se = seudo(getIm,P,'timeCourses',T,'nFrames',datInf.nFrames);                 % Load with NO classification
else
    se = seudo(getIm,P,'timeCourses',T,'nFrames',datInf.nFrames,'tcDefault',...
                                                     [datInf.Dir,datInf.classFile]); % Load with classification
end
% 'tcDefault','/home/adam/GITrepos/calcium_modeling/511510779_three_session_A_503109347_SEUDOclassification.mat')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Classify transients
%
% Instructions: 
%    Classify transients as "true" (blue) "false" (red) "mixed" (green) 
%    "unclassifiable" (gray) and "artifact" (brown). Definitions of these 
%    are:
%
%    True:           Only the cell that the profile outlines (even if it's 
%                    only half of the cell) is active
%    False:          A different source (potentially overlapping a lot with 
%                    the profile is active
%    Mixed:          Multiple sources, including the true source, are 
%                    active simultaneously
%    Unclassifiable: No discernable source, true or false, is active (the 
%                    image seems to be just noise)
%    Artifact:       Artifacts are profiles that clearly are capturing two
%                    or more different cells (i.e. half the profile lights
%                    up some times and the other half other times). E.g., 
%                    this happens with part of a cell and an apical 
%                    dendrite being lumped together. Artifacts are checked
%                    off in the bottom left hand corner of the GUI.
%
% Useful features:
%    - Can press 'o' to cycle through highlighting the overlap with shading, 
%      a boundary, or not at all
%    - Can press a number (0,1,2,...) to blur the image with a Gaussian 
%      filter of that width (in pixels). 0 = no blur. I've found 1 to be a
%      good value to use. This helps see weak transients over the noise
%    - The dynamic range can be changed to span only the pixels inside the 
%      profile. The dropdown menue controlling this is in the bottom-middle
%    
%  REMEBER TO SAVE YOUR CLASSIFICATIONS - AT THE END AND INTERMITTANTLY

se.classifyTransients

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

countMatrix = getClassCounts(se);

stats.All    = sum(sum(countMatrix));
stats.Nneur  = size(countMatrix,1);
stats.NArt   = sum(sum(countMatrix(:,1:3),2)==0);
stats.PerArt = 100*stats.NArt/stats.Nneur;
stats.TruRaw = sum(countMatrix(:,1));
stats.FalRaw = sum(countMatrix(:,2));
stats.MixRaw = sum(countMatrix(:,3));
stats.TruPer = 100*stats.TruRaw/sum(sum(countMatrix));
stats.FalPer = 100*stats.FalRaw/sum(sum(countMatrix));
stats.MixPer = 100*stats.MixRaw/sum(sum(countMatrix));
stats.ArtRaw = sum(countMatrix(:,end).*(sum(countMatrix(:,1:3),2)==0),1);
stats.ArtPer = 100*stats.ArtRaw/sum(sum(countMatrix));
stats.UncRaw = sum(countMatrix(:,end).*(sum(countMatrix(:,1:3),2)~=0),1);
stats.UncPer = 100*stats.UncRaw/sum(sum(countMatrix));

%%

histogram(countMatrix(:,3)./(sum(countMatrix,2).*(sum(countMatrix(:,1:3),2)>0)),20)
box off
% title('Experiment ID: 503109347')
title('Experiment ID: 506823562')
xlabel('Ratio of false transients')
ylabel('Number of cells')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 

% cellsSubset = randsample(293,30);
load('/mnt/bucket/adam/SEUDOdata/ALLEN_INST/506823562_subsetToClass.mat','cellsSubset')
mA = load('/home/adam/GITrepos/calcium_modeling/506823562_10pctClassAdam.mat', 'transientInfo');
mE = load([datInf.Dir,'adam_seudo_classification2.mat'], 'transientInfo');

seA.tcDefault.transientInfo = mA.transientInfo(cellsSubset);
seE.tcDefault.transientInfo = mE.transientInfo(cellsSubset);

cmA = getClassCounts(seA);
cmE = getClassCounts(seE);

clear mA mE

%%

tcCor     = se.cleanTimeCoursesWithClassification('default');
tcLSQ     = se.tcDefault.tc;
cov.tcCor = tcCor.'*tcCor - mean(tcCor,1).'*mean(tcCor,1);
cov.tcLSQ = tcLSQ.'*tcLSQ - mean(tcLSQ,1).'*mean(tcLSQ,1);
corr.tcCor = diag(sqrt(1./diag(cov.tcCor)))*cov.tcCor*diag(sqrt(1./diag(cov.tcCor)));
corr.tcLSQ = diag(sqrt(1./diag(cov.tcLSQ)))*cov.tcLSQ*diag(sqrt(1./diag(cov.tcLSQ)));

%%

histogram(mean((corr.tcCor - corr.tcLSQ)./corr.tcLSQ))

sum(sum(((corr.tcCor - corr.tcLSQ)./corr.tcLSQ)>0.1))./(numel(corr.tcLSQ)-size(corr.tcLSQ,1))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
