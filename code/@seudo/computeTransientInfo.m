function computeTransientInfo(se, whichStruct, varargin)
% seudoObject.computeTransientInfo(whichStruct,...)
%
% identify transients and compute summary information about each
%   (this function is called automatically when needed)
%
%
%
% EXAMPLE USAGE
% 
%   se.computeTransientInfo('default')
%   se.computeTransientInfo('contam','winRadius',14)
%   se.computeTransientInfo({'seudo',3},'threshFcn',@(x)5*std(x))
%
%
%
% INPUTS
%
%   whichStruct - specifies which time courses to compute transient info for
%                   e.g. 'default', 'contam', {'seudo',3}
%                   
%
% PARAMETERS
%
%
%   transientFrames - T x C logical matrix indicating when each cell exhibited a transient.
%                       if supplied, several other parameters are ignored (threshFcn, blurRadius, minDuration)
%
%   threshFcn - function applied to each cell time course to identify time points
%                   containing a significant increase in fluorescence. in the binarized
%                   time course, each connected component is considered a candidate transient.
%   blurRadius - candidate time courses are extended by length blurRadius in both directions
%                   before their onset and offset is computed. this serves to join components
%                   that are separated by small gaps.
%   minDuration - minimum duration (in frames) for a candidate transient to be included
%
%   tPre,tPost - # frames before transient onset and after transient offset to be used in analysis
%   winRadius - the spatial analysis is expanded in all directions by winRadius pixels 
%   useCOM - if true, winRadius is applied from the pixel at the profile's center of mass
%            if false, winRadius is applied from the minimum rectangle that includes
%                   all non-zero pixels in the profile
%
% 
%
% OUTPUT
%
%  	results are saved in a struct array, each element corresponding to one cell, with fields:
%   
%       shapes - Y x X x T matrix, each frame corresponds to one transient
%               each pixel is the weight of the transient time course when regressed against the pixel time course
%
%       corrWithProfile - T-length vector, correlation between transient shape and profile shape
%
%       shapesMean - same size as shapes
%               each pixel is the mean pixel time course during the transient
%
%       shapesCorr - same size as shapes
%               each pixel is the correlation between the transient time course and the pixel time course
%
%       times - T x 2 matrix, each row is the start and end frame of the transient
%
%       window - 4-length vector of the form [xStart xEnd yStart yEnd]
%                       specifies the spatial window used for fields shapes, shapesMean, and shapesCorr
%                       
%       classification - T-length vector indicating the classification of each transient
%
%                   1 = true
%                   0 = false
%                 0.5 = mixed
%                 nan = unclassified 
%
%               the classification is initialized to be all nans. if a classification already exists
%               when computeTransientInfo is called, the code attempts to preserve the classification for
%               as many transients as possible.
%       
%       params - struct, copy of the parameters used to compute transient info
%
%
%   
%       note: 
%
%           T = # transients
%           Y = # vertical pixels
%           X = # horizontal pixels
%
% 2018 Jeff Gauthier & Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETERS


% parse inputs
p = inputParser;
p.addParameter('transientFrames' , []); 
p.addParameter('threshFcn'         , @(x)max(seudo.robustSTD(x)*5,prctile(x,95))); 
p.addParameter('blurRadius'        , 2); 
p.addParameter('winRadius'         , 10);  
p.addParameter('useCOM'            , false);  
p.addParameter('minDuration'       , 6);
p.addParameter('tPre'              , 5); 
p.addParameter('tPost'             , 10); 
parse(p,varargin{:});
params = p.Results;
clear p


transientFrames   = params.transientFrames;
threshFcn         = params.threshFcn;
blurRadius        = params.blurRadius;
winRadius         = params.winRadius;
useCOM            = params.useCOM;
minDuration       = params.minDuration;
tPre              = params.tPre;
tPost             = params.tPost;

% put seudo object data into convenient variables
structSpec  = pickTcStruct(se,whichStruct);
profiles    = se.profiles;
timeCourses = se.(structSpec{1})(structSpec{2}).tc;
nCells      = size(profiles,3);

% ensure consistent sizes

assert(size(timeCourses,2) == size(profiles,3))
if ~isempty(transientFrames)
    assert(isequal(size(timeCourses),size(transientFrames)))
end

if ~isequal(size(profiles(:,:,1)),[se.movY se.movX])
    error('size of ROIs image frame (%d x %d) not equal to size of movie image frame (%d x %d)',size(profiles,1),size(profiles,2),se.movY,se.movX)
end

assert(size(timeCourses,1) == se.movF)

% DIVIDE INTO DISTINCT TRANSIENTS

% compute thresholded time courses, if not provided by user
if isempty(transientFrames)
    transientFrames = seudo.identifyTransients(timeCourses,...
        'threshFcn',threshFcn,'blurRadius',blurRadius,'minDuration',minDuration);
end

% IDENTIFY TRANSIENT START/STOP FRAMES AND SPATIAL WINDOW

transStartStopMatrix = []; % columns: [ (cell #) (transient start frame) (transient stop frame) ]
roiWindows = nan(4,nCells);


for cc=1:nCells
    
    % segment thresholded time course
    tcSeg = bwlabel(transientFrames(:,cc));
    
    % get # transients
    nTrans = max(tcSeg);
    
    % for each transient, note start/stop frame
    matThisCell = nan(nTrans,3);
    matThisCell(:,1) = cc;
    for tt=1:nTrans
        tBinary = tcSeg==tt;
        startFrame = max(1,find(tBinary,1,'first')-tPre);
        stopFrame = min(se.movF,find(tBinary,1,'last')+tPost);
        matThisCell(tt,2:3) = [startFrame stopFrame];
    end
    
    % add to accumulating list 
    transStartStopMatrix = cat(1, transStartStopMatrix, matThisCell);
    
    
    
    % identify spatial window
    
    % get outer bounds of area to which to add padding
    if useCOM % just the COM
        theCom = seudo.computeRoiCOMs(profiles(:,:,cc));
        ob = round(theCom([2 2 1 1]));
    else % outer bound of non-zero pixels
        [~,ob] = seudo.computeRoiCOMs(profiles(:,:,cc));
    end
    
    % add padding
    roiWindows(:,cc) = [max(1,ob(1)-winRadius) min(se.movY,ob(2)+winRadius)  max(1,ob(3)-winRadius) min(se.movX,ob(4)+winRadius)]';
    
end


% initialize final storage variable
transientInfo = struct('shapes',cell(nCells,1),'shapesCorr',cell(nCells,1),'shapesMean',cell(nCells,1));

for cc=1:nCells
    % load basic info
    transientInfo(cc).times = transStartStopMatrix(transStartStopMatrix(:,1)==cc,2:3);
    transientInfo(cc).window = roiWindows(:,cc)';
            
    % initialize classification
    transientInfo(cc).classification = nan(size(transientInfo(cc).times,1),1);
end


% procedure:
%
% load new frame
% for each cell with a transient in this frame
%      put frame subset into the current mini-movie (create if needed)
%      if a mini-movie is now complete
%           do active shape analysis on it
%           delete mini-movie
%



% initialize storage of mini-movies, each centered on the cell profile
transientMovies = cell(size(transStartStopMatrix,1),1);


% create lookup matrix for which transients take place during which frames

% transLookupMat = logical(sparse(size(transStartStopMatrix,1),se.movF));
transLookupCell = cell(se.movF,1);
for tt=1:size(transStartStopMatrix,1)
    for ff = transStartStopMatrix(tt,2):transStartStopMatrix(tt,3)
        transLookupCell{ff} = cat(2, transLookupCell{ff}, tt);
    end
end



% identify which frames have a transient
transLookupMat = false(size(timeCourses));
for tt=1:size(transStartStopMatrix,1)
    transLookupMat(transStartStopMatrix(tt,2):transStartStopMatrix(tt,3),transStartStopMatrix(tt,1)) = true;
end
transFrames = find(sum(transLookupMat,2)>0)';



% confirm result has not changed from slow version
% if 0 || 0
%     transLookupMat_slow = sparse(false(size(transStartStopMatrix,1),se.movF));
%     for tt=1:size(transStartStopMatrix,1)
%         transLookupMat_slow(tt,transStartStopMatrix(tt,2):transStartStopMatrix(tt,3)) = true;
%     end
%     assert(isequal(find(sum(transLookupMat_slow,1)),transFrames))
% end



% for every frame that has a transient in it...
for ff = transFrames
    
    
    % load the frame
    try
    theFrame = se.getFrame(ff);
    catch ME
        fprintf('Unable to load frame %d\n',ff)
        rethrow(ME)
    end
    
    % identify which transients are ongoing during this frame
    transOngoing = transLookupCell{ff};
    
    % confirm result has not changed from slow version
%     if 0 || 0
%         transOngoing_slow = find(transLookupMat_slow(:,ff))';
%         assert(isequal(transOngoing_slow,transOngoing))
%     end
    
    
    % accumulate list of transients ending at this frame
    transStopping = [];
    
    % add this frame to mini-movie of each ongoing transient
    for tt = transOngoing
        
        % get window
        cutWin = roiWindows(:,transStartStopMatrix(tt,1));
        
        % note which was the first frame of this transient
        startFrame = transStartStopMatrix(tt,2);
        
        % initialize, if first time
        if startFrame == ff
            transientMovies{tt} = nan(cutWin(2)-cutWin(1)+1,cutWin(4)-cutWin(3)+1,transStartStopMatrix(tt,3)-transStartStopMatrix(tt,2)+1);
        end
        
        % insert subset of ROI into this mini-movie
        transientMovies{tt}(:,:,ff-startFrame+1) = theFrame(cutWin(1):cutWin(2),cutWin(3):cutWin(4));
        
        % note if ending
        if transStartStopMatrix(tt,3) == ff
            transStopping = cat(2, transStopping, tt);
        end
    end
    
    
    % compute active shape for any transients that just ended
    if ~isempty(transStopping)
        for tt = transStopping
            
            % get transient time course
            tc = timeCourses(transStartStopMatrix(tt,2):transStartStopMatrix(tt,3),transStartStopMatrix(tt,1));
            
            % use mini-movie to compute active shape
            results = computeOneActiveShape(tc,transientMovies{tt});
            
            % delete mini-movie
            transientMovies{tt} = [];
            
            % store results in struct
            cc = transStartStopMatrix(tt,1);
            transientInfo(cc).shapes     = cat(3,transientInfo(cc).shapes,results.shape);
            transientInfo(cc).shapesCorr = cat(3,transientInfo(cc).shapesCorr,results.shapeCorr);
            transientInfo(cc).shapesMean = cat(3,transientInfo(cc).shapesMean,results.shapeMean);
        end
    end
    
end






% if there is an existing classification, try to resolve it with the new set of transients
if isfield(se.(structSpec{1})(structSpec{2}),'transientInfo') && ...
        isfield(se.(structSpec{1})(structSpec{2}).transientInfo,'classification')
    
    % ensure the number of cells is the same in the old and new transient info structs
    assert(length(se.(structSpec{1})(structSpec{2}).transientInfo) == length(transientInfo))
    
    for cc = 1:nCells
        
        % copy the artifact classification
        transientInfo(cc).isArtifact = se.(structSpec{1})(structSpec{2}).transientInfo(cc).isArtifact;
        
        % get the list of old and new transient times
        oldTimes = se.(structSpec{1})(structSpec{2}).transientInfo(cc).times;
        newTimes = transientInfo(cc).times;
        
        % if identical, just copy the classification
        if isequal(oldTimes,newTimes)
            transientInfo(cc).classification = se.(structSpec{1})(structSpec{2}).transientInfo(cc).classification;
        else
            
            % if not identical, try to resolve differences
            
            % make matrix of old transient times
            tcOld = zeros(se.movF,nCells);
            for tt = 1:size(oldTimes,1), tcOld(oldTimes(tt,1):oldTimes(tt,2),tt) = tt; end
            
            % for each new transient...
            for tt = 1:size(newTimes,1)
                
                % identify which old transients overlap this new transient
                whichOldOverlap = unique(tcOld(newTimes(tt,1):newTimes(tt,2),:));
                whichOldOverlap = setdiff(whichOldOverlap,0);
                
                
                switch length(whichOldOverlap)
                
                    case 0 % if none overlap, leave unclassified (default value)
                        
                    case 1  % if just one old transient overlaps, copy its classification
                        transientInfo(cc).classification(tt) = ...
                            se.(structSpec{1})(structSpec{2}).transientInfo(cc).classification(whichOldOverlap);
                        
                    otherwise % if multiple old transients overlap, combine them
                        oldVals = num2cell(se.(structSpec{1})(structSpec{2}).transientInfo(cc).classification(whichOldOverlap));
                        
                        % if they're all the same, keep that value
                        if isequaln(oldVals{:})
                            transientInfo(cc).classification(tt) = oldVals{1};
                        else
                            % otherwise set to unclassified
                            transientInfo(cc).classification(tt) = seudo.valUnc;
                        end
                end
                
            end
        end
    end
    
    
end


% if there is an existing artifact classification, copy it
if isfield(se.(structSpec{1})(structSpec{2}),'transientInfo') && ...
        isfield(se.(structSpec{1})(structSpec{2}).transientInfo,'isArtifact')
    
    % copy the artifact classification
    for cc = 1:nCells
        transientInfo(cc).isArtifact = se.(structSpec{1})(structSpec{2}).transientInfo(cc).isArtifact;
        
    end
else
    
    % otherwise initialize the field for artifact classification
    for cc = 1:nCells
        transientInfo(cc).isArtifact = false;
    end
end





% get correlation between cell profile and each transient
for cc = 1:length(transientInfo)
    
    % get profile
    prof = se.profiles(transientInfo(cc).window(1):transientInfo(cc).window(2),transientInfo(cc).window(3):transientInfo(cc).window(4),cc);
    
    % compute correlation of each shape
    if ~isempty(transientInfo(cc).shapes)
        corrVals = correlationVectorMatrix(prof(:),reshape(transientInfo(cc).shapes,[],size(transientInfo(cc).shapes,3)))';
    else
        corrVals = [];
    end
    
    % store
    transientInfo(cc).corrWithProfile = corrVals;
    
    % add a copy of the parameters
    transientInfo(cc).params = params;
end






% store result in seudo object
se.(structSpec{1})(structSpec{2}).transientInfo = transientInfo;

end







function results = computeOneActiveShape(tc,mov)



% reshape movie matrix
wY = size(mov,1);
wX = size(mov,2);
movReshape = reshape(mov,wY*wX,length(tc)); % pix X time

if all(isnan(tc))
    
    % if there is no time course, fill in with NaNs
    fitToTC = nan(wY,wX);
    corrWithTC = nan(wY,wX);
    meanImage = nan(wY,wX);
    
else
    % identify NaNs in the time course (which shouldn't be there anyway, but here we are)
    safePoints = ~isnan(tc);
    
    % compute how well each pixel matches the time course
    fitToTC = reshape(movReshape(:,safePoints) / tc(safePoints)',wY,wX);
    
    % compute correlation of each pixel with extracted time course
    corrWithTC = reshape(correlationVectorMatrix(tc,movReshape'),wY,wX);
    
    % compute mean image
    meanImage = mean(mov,3);
end


results = struct;

results.shape = fitToTC;
results.shapeCorr = corrWithTC;
results.shapeMean = meanImage;
end
