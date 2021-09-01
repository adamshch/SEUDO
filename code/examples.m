% This file shows example usage of most methods. It should not be executed in full,
% since several equivalent variants of each command are shown.
%
% For a demonstration using a provided dataset, see demo.m
% For list of properties in the seudo object class, see documentation in seudo.m
% For full list of methods, see seudo.m and files in folder @seudo.
%
%

%#ok<*NASGU>




% LOAD DATA INTO SEUDO OBJECT


% ----> if your whole movie fits into memory

% provide
%   1) dF/F movie as a matrix M, size [ Y x X x (# frames) ]
%   2) source profiles as a matrix P, [ Y x X x (# sources) ]
se = seudo(M,P);



% ----> if your whole movie does NOT fit into memory

% provide 
%   1) dF/F movie in one of the formats below
%   2) source profiles as a matrix P, [ Y x X x (# sources) ]

% specify dF/F movie as tiff files
se = seudo({'/data/ims001.tif','/data/ims002.tif','/data/ims003.tif'},P);

% specify movie as tiff files that contain absolute flourescence, and also provide a baseline image 
% and zero level. the baseline image will be subtracted and divided out to yield dF/F)
se = seudo({'/data/ims001.tif','/data/ims002.tif','/data/ims003.tif'},P,'baseline');

% specify movie as user function that takes frame number as its only input
se = seudo(@getMovieFrame,P);

% specify movie as user function that takes frame number in addiiton to other inputs
se = seudo( @(whichFrame) loadMovieFrame('mouse_006','June 27',whichFrame) ,P);




% OPTIONAL THINGS TO LOAD

% provide time courses
se = seudo(M,P,'timeCourses',T);

% load the whole movie into a large matrix, then use it in different SEUDO objects
M = bigMovieMatrix;
se_001 = seudo( @(f)M(:,:,f), P_001);
se_002 = seudo( @(f)M(:,:,f), P_002);

% load without specifying movie, but supplying time course struct elsewhere
se = seudo([],P,'tcDefault','/analysis/savedClassification.mat');



% OPTIONAL provide transient definitions
se = seudo(M,P,'timeCourses',T);



% PERFORM SEUDO

% show suggested parameters (this is also automatically done when loading the dataset)
se.suggestParameters

% default parameters (probably won't work!)
se.estimateTimeCoursesWithSEUDO

% manually specified parameters
se.estimateTimeCoursesWithSEUDO('blobRadius',2.5,'lambdaBlobs',0.01)








% OPTIONAL: IDETNIFY OPTIMAL PARAMETERS, QUICK VERSION

% choose a few cells and transients
se.classifyTransients('default')

% put selections in a matrix: each row is cell number then transient number, e.g.
pickedCellTransients = [ 4 3 true; 113 6 false; 134 9 false];

% perform SEUDO in real time on these
%   launches a GUI with text field to enter parameters
%   hit return, and it performs SEUDO on these transients
se.quickSEUDO('default',pickedCellTransients)

% apply less spatial padding around each profile, and include more frames around each transient
se.quickSEUDO('default',pickedCellTransients,'padSpace',0,'padTime',20)






% OPTIONAL: IDENTIFY OPTIMAL PARAMETERS, LONG VERSION

% apply the SEUDO aogorithm using several different parameter sets
%   each profile is considered separately, same as high contam time courses
result = se.paramSearch(pickedCellTransients,'blobRadius',[1 1.5 2],'p',[0.0001 0.01]);

% see results to find best match with manually classified transients
seudoReviewParamSearch(result)

% store best parameters in a struct
bestParams = result.paramSets{7};

% apply SEUDO algorithm using these parameters
se.estimateTimeCoursesWithSEUDO(bestParams{:})







% OPTIONAL: EVALUATE SEUDO STATISTICALLY 


% run for several different parameter sets
%   all versions of SEUDO analysis are loaded into the seudo object
se.parallelSEUDO([], result.paramSets{6}{:})
se.parallelSEUDO([], result.paramSets{7}{:})
se.parallelSEUDO([], result.paramSets{8}{:})

% plot ROC curve
se.computeSeudoROC







% OPTIONAL: COMPARE SEUDO TO MAXIMALLY CONTAMINATED TIME COURSES

% compute contaminated time courses (might take several minutes or hours)
% this ignores all other found cells to create many false transients
se.computeLSQtimeCoursesHighContam

% manually classify contaminated transients
se.classifyTransients('contam')

% perform SEUDO using the same high contamination constraint
se.estimateTimeCoursesWithSEUDO(bestParams{:},'contam',true)

% compare the results with and without seudo
se.compareTimeCourses('contam','seudo')




