%% LOAD ALLEN DATA INTO SEUDO OBJECT


%% CHOOSE DATASET


dataDir = '/Volumes/tank/jlgauthi/SEUDO/ALLEN_INST/';

switch 1
    case 1
        fnIms = '511510779_three_session_A_503109347.h5';
        fnT = 'allendata_dff_503109347.h5';
    case 2
        fnIms = '511510779_three_session_B_503019786.h5';
        fnT = 'allendata_dff_503019786.h5';
    case 3
        fnIms = '511510779_three_session_C_502634578.h5';
        fnT = 'allendata_dff_502634578.h5';
    case 0
        
%         '511510848_three_session_A_506823562.h5'
%         '511510848_three_session_B_507691834.h5'
%         '511510848_three_session_C_504809131.h5'
%         
%         '536323956_three_session_A_540684467.h5'
%         '536323956_three_session_B_541048140.h5'
%         '536323956_three_session_C_539515366.h5'
%         
%         '561472631_three_session_A_561472633.h5'
%         '561472631_three_session_B_561528269.h5'
%         '561472631_three_session_C2_562003583.h5'
%         
%         '562536151_three_session_A_562536153.h5'
%         '562536151_three_session_B_567446262.h5'
%         '562536151_three_session_C2_563027941.h5'
%         
%         '564425775_three_session_A_564425777.h5'
%         
%         
%         'allendata_dff_504809131.h5'
%         'allendata_dff_506823562.h5'
%         'allendata_dff_507691834.h5'
%         'allendata_dff_539515366.h5'
%         'allendata_dff_540684467.h5'
%         'allendata_dff_541048140.h5'
%         'allendata_dff_561472633.h5'
%         'allendata_dff_561528269.h5'
%         'allendata_dff_562003583.h5'
%         'allendata_dff_562536153.h5'
%         'allendata_dff_563027941.h5'
%         'allendata_dff_564425777.h5'
%         'allendata_dff_567446262.h'
end
fnP = strrep(fnT,'dff','rois');

%% GET SIZES

H = h5info([dataDir fnIms],'/data');
nY = H.Dataspace.Size(1);
nX = H.Dataspace.Size(2);
nFrames = H.Dataspace.Size(3);
%nFrames = 10000;

H = h5info([dataDir fnT],'/dff');
nProfiles = H.Dataspace.Size(2);

clear H

%% GET MEAN IMAGE

dataForMean = h5read([dataDir fnIms],'/data',[1 1 1],[nY nX 100],[1 1 floor(nFrames / 100)]);
meanIm = double(mean(dataForMean,3));

%% GET TIME COURSES & ROIs

T = h5read([dataDir fnT],'/dff',[1 1],[nFrames nProfiles],[1 1]);
P = double(h5read([dataDir fnP],'/X',[1 1 1],[nY nX nProfiles],[1 1 1]));


%% SET UP FUNCTION

getIm = @(f)double(h5read([dataDir fnIms],'/data',[1 1 f],[512 512 1],[1 1 1]))-meanIm;

se = seudo(getIm,P,'timeCourses',T,'nFrames',nFrames,'name',strrep(fnIms,'.h5',''));


%% LOAD EXISTING CLASSIFICATION

se = seudo(getIm,P,'timeCourses',T,'nFrames',nFrames,strrep(fnIms,'.h5',''),'tcDefault',...
    '/Volumes/tank/jlgauthi/SEUDO/Allen_classification/511510779_three_session_A_503109347_JG_inc_class.mat');

