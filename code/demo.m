% execute this file in sequential blocks to see a demonstration of how to use SEUDO with a provided dataset



%%



%%%%%%%%%%%%%%%%%%%%%
%                   %
%     LOAD DATA     %
%                   %
%%%%%%%%%%%%%%%%%%%%%




%% LOAD .mat FILE

% a demo dataset can be downloaded here:
%    https://drive.google.com/file/d/1N-ESEm1eU2G6wRVd9LdZ8bOLqmSUqTgg

% you can either drag demoData.mat onto the matlab command window, or execute
% this block after setting the current path to the file location

load('demoData1.mat')


%% PACKAGE INTO SEUDO OBJECT AND PROVIDE CUSTOM TRANSIENT DEFINITIONS

% M - movie data
% P - source profiles
% T - source time courses
% TF - logical matrix specifying which sources had transients in which frames

se = seudo(M,P,'timeCourses',T,'name','demo dataset');
se.computeTransientInfo('default','transientFrames',TF)




%%



%%%%%%%%%%%%%%%%%%%%%
%                   %
%   EVALUATE DATA   %
%                   %
%%%%%%%%%%%%%%%%%%%%%




%% OPTIONAL: LOAD OUR CLASSIFICATION INTO THE SEUDO OBJECT
% you can load our classification if you don't want to manually classify transients

se = seudo(M,P,'timeCourses',T,'name','demo dataset','tcDefault',tDefaultDemo);



%% LAUNCH TRANSIENT INSPECTION/CLASSIFICATION GUI
% learn more about the interface in the documentation.pdf, section
% "Look through transients to classify them manually"

se.classifyTransients



%% FIND GOOD PARAMETERS FOR AUTOMATIC CLASSIFICATION
% this uses the current classification as "ground truth", so for best results be sure
% that many transients are classified

% auto classification is based on a correlation threshold (see manuscript for details),
% and here several possible thresholds are compared

% takes ~20 seconds to compute
se.optimizeAutoClassification('threshValues',-.5:.1:1)



%% MAKE FALCON PLOT
% apply automatic classification to estimate what fraction of transients are true for each source,
% colors indicate the human classification for comparison

se.plotFalcon([],[],'corrThresh',0.5)




%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                          %
%   RUN SEUDO TO REMOVE FALSE TRANSIENTS   %
%                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%% GET A BALLPARK ESTIMATE FOR PARAMETERS FOR SEUDO
% currently this only generates a suggestion for one parameter, 'sigma2'

se.suggestParameters



%% RAPIDLY TEST SEUDO PARAMETERS ON A FEW TRANSIENTS
% this GUI helps you identify parameters for optimal performance, i.e. preserving
% true transients and removing false transients


% specify a small number of transients (here one true and one false) with format:
%    [ sourceID  transientID  transientClass ; ... ]
selectedTransients = [50 61 true; 50 74 false];

% launch GUI (will take 10-20 seconds to compute results before displaying)
se.quickSEUDO(selectedTransients)

% try changing parameter 'lambdaBlob' to see how the results change



%%  COMPARE SEUDO PERFORMANCE OVER A RANGE OF PARAMETERS
% after identifying parameters near the optimal regime, you can try out a wide range of
% parameter combinations, then view the results in a GUI


% inspect the results
seudoReviewParamSearch(psResults)


% to recompute the result, run the code in this sub-section
if 0 
    %% it takes about an hour to run on our system
    cellSpec = [12 42 false; 12 7 true; 50 61 true; 50 74 false];
    psResults = se.paramSearch(cellSpec,'factorial',...
        'sigma2',[0.02 0.2 2],'lambdaBlob',[0.01 0.1 1 10],'dsTime',[1 3],'blobRadius',[2 3.5 5]);
end



%% OPTIMIZE PARAMETERS WITH BAYESIAN ALGORITHM
% you can also find optimal parameters using BayesOpt

% requires Matlab 2016b or later

cellSpec        = [50 61 true; 50 74 false];
bayesOptResults = se.bayesianParameterOptimization(cellSpec);






%%




%%%%%%%%%%%%%%%%%%%%%
%                   %
%     RUN SEUDO     %
%                   %
%%%%%%%%%%%%%%%%%%%%%





%% LOAD .mat FILE


% the second demo file can be downloaded here:
%    https://drive.google.com/file/d/1u6aSlN9vtOY9DZa0PJVFmc9iqKQaOFUN

% NOTE: this variable was saved in a sepaerate file because it occupies 8 GB of memory

load('demoData2.mat')


%% COMPUTE TIME COURSES WITH SEUDO
% SEUDO is computed for each source separately, and only in frames with a transient. Computation takes about
% 2-10 minutes per source.


% select parameters
seudoParams = {'dsTime',3,'sigma2',0.0020,'lambdaBlob',10,'blobRadius',3,'padSpace',5,'saveBlobTimeCourse',1};
%seudoParams = {'dsTime',3,'sigma2',0.0020,'lambdaBlob',300,'blobRadius',3,'padSpace',5,'saveBlobTimeCourse',1};

% select cells
whichCells = 1:10;

% run SEUDO
se.parallelSEUDO(whichCells, seudoParams{:})





%% COMPARE RESULTS FROM SEVERAL PARAMETER SETS

% run for a second set of parameters
seudoParams = {'dsTime',3,'sigma2',0.0020,'lambdaBlob',1,'blobRadius',3,'padSpace',5,'saveBlobTimeCourse',1};
se.parallelSEUDO(whichCells, seudoParams{:})

% show results
se.computeSeudoROC


% optional: save ROC numbers to a variable
%theResults = se.computeSeudoROC;

