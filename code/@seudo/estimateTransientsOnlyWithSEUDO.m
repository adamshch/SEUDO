function se = estimateTransientsOnlyWithSEUDO(se,varargin)
% seudoObject.estimateTimeCoursesWithSEUDO(...)
%
% estimate time courses using SEUDO algorithm
%
%
% EXAMPLE:
%
%   seudoObject.estimateTimeCoursesWithSEUDO('p',0.001','lambdaBlob',1)
%
%
%
% PARAMETERS:
%
%   p - probability of false transient
%   sigma2 - per-pixel noise variance
%   lambdaBlob - sparsity multiplicative value for blobs
%   blobRadius - radius (in pixels) of a single Gaussian blob
%
%   dsTime - # frames to average. SEUDO analysis is applied to the average of the most recent <dsTime> frames 
%   lambdaProf - sparsity multiplicative value for found cell profiles
%
%   minPixForInclusion - how many pixels must be in another found ROI to be included in each cell's SEUDO. Minimum value is 1.
%   padSpace - defines width around the cell's minimum rectangle in which pixels are used for SEUDO
%   useCOM - instead of minimum rectangle, use the cell's COM as the point around which the padSpace width applies
%
%   whichCells - list of which cells to perform analysis on. if empty (default), use all.
%                   unselected cells have a time course of NaNs.
%
%   verbose - logical, whether to print progress to the command window
%   showDefaults - if true, indicate default values without performing SEUDO
%                       if no output variables are set, just display the values
%                       otherwise, put a struct of default values into the output variable 
%   validateParams - if true, validate parameters without performing SEUDO
%                       if all parameters are valid, returns empty
%                       if any parameters are invalid, returns the first error
%
%   saveBlobTimeCourse - logical, whether to save the fit weight of every blob in every frame
%                           for very large movies, setting to false (the default) can save memory
%
%   saveOtherCellTimeCourses - logical, whether to save the fit weight of cells other than the cell of interest
%                           for very long movies with many cells, setting to false (the default) can save memory
%   
%   dyn_mod - string specifying which algorithm to use
%               'static' - single frame
%               'rwl1df' - 
%               'bpdndf' - 
%
%   parameters only affecting 'bpdndf'
%
%       kappa - 
%       kappaProf - 
%
%   parameters only affecting 'rwl1df' (re-weighted L1)
%
%       n_rw_iter - 
%       tau - 
%       eta - 
%       beta - 
%
%   
%
%
%
%
% OUTPUTS:
%
%   All results are stored in the seudo object in field 'tcSeudo'. This field is a struct array
%   in which each element is one set of results. The most recent results are stored in the last
%   struct, i.e. seudoObject.tcSeudo(end)
%
%   Each struct contains these fields:
%
%       tc - (# frames) x (# cells) matrix, each column the time course of one found cell
%       tcLSQ - (# frames) x (# cells) matrix, each column the least squares time course of one found cell
%       params - struct with a copy of the paramters used to run SEUDO,
%                   plus field 'version' indicating the software version 
%       extras - struct storing intermediate results in the computation, with these fields:
%
%               lsqCosts - (# frames) x (# cells) matrix, cost of least squares solution in each frame
%
%               tcBlobs - time courses of blobs. format depends on parameter saveBlobTimeCourse.
%                   if true: {# cells} (# frames) x (# pixels analyzed for each cell) matrix
%                           each column is the time course of one blob
%                 	if false:  (# frames) x 2 x (# cells) matrix
%                           first column is sum of all blobs, second column is sum of blobs within the ROI
%
%               whichPixels - 
%               whichProfiles - 
%               tcOtherCellsLSQ - 
%
%               ...
%   
%
%
% 2019 Adam Chrles & Jeff Gauthier


if ~exist('tfocs','file')
    error('TFOCS must be installed to use SEUDO, download from http://cvxr.com/tfocs/download/')
end
CODE_VERSION = 1;                                                          % Set code version

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input Parameter parsing

% Get Input Parameters
p = setDefaultInputs();
q = makeValidationParser();
parse(q,varargin{:});

if q.Results.validateParams                                                % just validate inputs, if desired
    se = [];                                                               % if inputs can be parsed successfully, return empty
    try parse(p,varargin{:}); catch err, se = err; end                     % otherwise return the error
    return
end

if q.Results.showDefaults                                                  % just show defaults, if desired
    parse(p);                                                              % get default values
    params         = p.Results;
    params.version = CODE_VERSION;
    if nargout > 0                                                         % if there are output variables...
        se = params;                                                       %  ...return struct
        return
    else                                                                   % otherwise, just display default values
        fprintf('\ndefault values for SEUDO:\n\n')
        disp(params)
        clear se
        return
    end
end

parse(p,varargin{:});                                                      % parse inputs
params         = p.Results;                                                % copy to struct that will be output
params.version = CODE_VERSION;                                             % set code version

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A bunch of initializations (variables & parameters)

if isempty(params.whichCells); whichCells = 1:size(se.profiles,3);
else;                          whichCells = params.whichCells(:)';
end

nFrames  = se.movF;                                                        % Number of time-frames
dsT      = params.dsTime;                                                  % temporal downsample
sigma2   = params.sigma2/dsT;                                              % noise variance
params.minPixForInclusion = max(params.minPixForInclusion,1);              % ensure at least 1
tc       = initTimeCourseCells(whichCells, nFrames);                       % initialize storage of final results
one_blob = makeOneBlob(params.blobRadius);                                 % make one blob

if params.saveBlobTimeCourse; tcBlobs = cell(length(whichCells),1);
else;                         tcBlobs = nan(nFrames,2,length(whichCells)); % save summed blob time course
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% initialize storage of cell-specific variables (1 cell index = 1 cell)

[md, whichProfiles, whichPixels, tc.tcBlobs] = collectMetaInfo(se, ...
                        whichCells, tc.tcBlobs, one_blob, sigma2, params); % Get a bunch of meta data containing the locations/pixels/parameters for SEUDO for each cell 
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get blob fit to each cell

blobFits = fitBlobsToProfiles(se, whichCells, params.padSpace, params.blobRadius);  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run per-frame analysis
 
updateFrames = createFramesToUpdateOn(se,[1 3 10 30 100 300 1000],params); % Create frames in which to update the user on progress
frameSlidingWindow = nan(se.movY*se.movX,dsT);                             % initialize holder of frames
startTime          = now;                                                  % Initialize start time
saveBlobsOpt = params.saveBlobTimeCourse;
saveOtherOpt = params.saveOtherCellTimeCourses;
dyn_mod       = params.dyn_mod;                                            % Temporary storage of the type of dynamical regularization to use
% tcBlobs = tc.tcBlobs;
tcTemp  = cell(length(whichCells),1);                                      % Initialize a place to store temporary results (to better slice variables for the parfor loop)

for ff=1:nFrames
    
    % get current frame
    currFrame = reshape(se.getFrame(ff),[],1);
    
    if ff <= dsT; frameSlidingWindow(:,ff) = currFrame;                    % the first few times, build up a matrix of the initial frames
    else; frameSlidingWindow = [frameSlidingWindow(:,2:end) currFrame];    % subequent times, delete the oldest frame and appent the current frame
    end
    thisFrameFull = nanmean(frameSlidingWindow,2);                         % average the last few frames
    if ff > 1; tcPrevSEUDO = tc.tcSEUDO(ff-1,:).';
    else;      tcPrevSEUDO = zeros(size(tc.tcSEUDO(ff,:).'));
    end
    for cc = 1:length(whichCells)                                       % Iteratively apply SEUDO to each cell
        if ff==1 
            tcTemp{cc}.fitBOB = []; 
            tcTemp{cc}.fitWeights = [];
        end
        tcTemp{cc} = runOneSEUDO(ff==1, md{cc}, thisFrameFull, ...
                            whichProfiles(:,cc), whichPixels(:,cc), ...
              tcPrevSEUDO,tcTemp{cc}.fitBOB,dyn_mod, params.optimizerOpt); % Run SEUDO on this frame for this cell
    end
    for cc = 1:length(whichCells) 
        % store fits in time course variables
        tc.tcSEUDO(ff,cc)          = tcTemp{cc}.tcSEUDO;
        tc.tcCellsWithBlobs(ff,cc) = tcTemp{cc}.tcCellsWithBlobs;
        tc.tcLSQ(ff,cc)            = tcTemp{cc}.tcLSQ;
                   
        if saveOtherOpt
            tc.tcOtherCellsLSQ_all{cc}(ff,:) = tcTemp{cc}.tcLSQ_thisFrame;
        end
        
        % save blobs
        if saveBlobsOpt
            tcBlobs{cc}(ff,:) = tcTemp{cc}.fitBOB;
        else
            tcBlobs(ff,1,cc) = sum(tcTemp{cc}.fitBOB);
            prof             = se.profiles(:,:,cc) > 0;
            tcBlobs(ff,2,cc) = sum(tcTemp{cc}.fitBOB(prof(whichPixels(:,cc))));
        end
        tc.lsqCosts(ff,cc) = tcTemp{cc}.lsqCosts;                          % store the Least-squares costs
        tc.bobCosts(ff,cc) = tcTemp{cc}.bobCosts;                          % store the blob costs
    end
    
    if ismember(ff,updateFrames)
        finishDur = (nFrames-ff)/ff*(now-startTime);
        fprintf('      completed %d of %d (%0.1f%%) in %0.1f min, estimated finish in %0.1f min (%s)\n',...
            ff, nFrames, 100*ff/nFrames,(now-startTime)*24*60,finishDur*24*60,datestr(now+finishDur))
        %whos
    end
    
end

% compile extras struct
extras                  = struct;
extras.lsqCosts         = tc.lsqCosts;
extras.blobCosts        = tc.bobCosts;
extras.tcCellsWithBlobs = tc.tcCellsWithBlobs;
extras.tcBlobs          = tcBlobs;
extras.one_blob         = one_blob;
extras.blobFits         = blobFits;
extras.whichPixels      = whichPixels;
extras.whichProfiles    = whichProfiles;
extras.tcOtherCellsLSQ  = tc.tcOtherCellsLSQ_all;


% STORE RESULTS IN SEUDO OBJECT
idx = length(se.tcSeudo) + 1;                                              % identify index to store in
se.tcSeudo(idx).tc     = tc.tcSEUDO;                                       % SEUDO time courses 
se.tcSeudo(idx).tcLSQ  = tc.tcLSQ;                                         % LSQ time courses 
se.tcSeudo(idx).params = params;                                           % parameters
se.tcSeudo(idx).extras = extras;                                           % intermediate computations for debugging

if params.verbose
    fprintf('   completed in %0.1f min (%s), results stored in field tcSeudo(%d)\n',...
        (now-startTime)*24*60,datestr(now),idx)
end

end

%% END OF MAIN FUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Extra functions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function p = setDefaultInputs()

p = inputParser;

p.addParameter('p'                  , 0.00001, @(x)x>=0&&x<=1);   
p.addParameter('sigma2'             , .01);   
p.addParameter('lambdaBlob'         , 20); 
p.addParameter('blobRadius'         , 1.2);

p.addParameter('dsTime'             , 1);
p.addParameter('lambdaProf'         , 0);
p.addParameter('minPixForInclusion' , 1);
p.addParameter('padSpace'           , 10);
p.addParameter('useCOM'             , false);

p.addParameter('verbose'        , true);
p.addParameter('validateParams' , false);
p.addParameter('showDefaults'   , false);

p.addParameter('dyn_mod'   , 'static');
p.addParameter('kappa'     , 0.05);
p.addParameter('kappaProf' , 0.001);
p.addParameter('n_rw_iter' , 3);
p.addParameter('tau'       , 0.1);
p.addParameter('eta'       , 0.01);
p.addParameter('beta'      , 0.01);

p.addParameter('whichCells'               , []);
p.addParameter('saveBlobTimeCourse'       , false);
p.addParameter('saveOtherCellTimeCourses' , true);

p.addParameter('optimizerOpt','TFOCS');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function q = makeValidationParser()

q               = inputParser;
q.KeepUnmatched = true;
q.addParameter('validateParams', false, @(x)isscalar(x)&&(islogical(x) || ismember(x,[0 1]))); % just validate inputs, if desired
q.addParameter('showDefaults' , false, @(x)isscalar(x)&&(islogical(x) || ismember(x,[0 1]))); % just show defaults, if desired
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create a single blob

function one_blob = makeOneBlob(blobRad)

% parameters
cropRad     = ceil(blobRad*2.5+eps);
clipHeight  = 0.01;

% make one blob and matrix of blobs
one_blob = fspecial('gauss',(cropRad*2+1)*[1 1],blobRad);
one_blob = one_blob .* double(one_blob > clipHeight * max(one_blob(:)));
one_blob = one_blob/sqrt(sum(one_blob(:).^2));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Function to fit each profile with a set of blobs

function blobFits = fitBlobsToProfiles(se, whichCells, padSpace, blobRad)

one_blob = makeOneBlob(blobRad);                                           % Make a single blob
cropRad  = ceil(blobRad*2.5+eps);                                          % Calculate the cropping radius
blobFits = nan(se.movX*se.movY,length(whichCells));                        % Initialize the matrix of blob fits

for cc = 1:length(whichCells)
    thisCell = whichCells(cc);                                             % get cell profile
    prof     = se.profiles(:,:,thisCell);                                  % get cell profile
    [~,ob]   = seudo.computeRoiCOMs(prof);                                 % identify region around profile
    ob       = [max(1,ob(1)-padSpace) min(se.movY,ob(2)+padSpace) ...
                       max(1,ob(3)-padSpace) min(se.movX,ob(4)+padSpace)]; % add padding equal to bblobb radius
   
    thesePix = false(se.movY,se.movX);                                     % identify which subset of the original this is
    thesePix(ob(1):ob(2),ob(3):ob(4)) = true;

    profSmall = prof(ob(1):ob(2),ob(3):ob(4));                             % get profile in this region
    blobMat_  = convmtx2(one_blob,size(profSmall));                        % construct matrix with blobs at all locations

    keepRegion = false(cropRad*2+size(profSmall));                         % clip just to the region where the whole blob is present
    keepRegion(cropRad+1:end-cropRad,cropRad+1:end-cropRad) = true;        % 
    
    blobMat  = blobMat_(keepRegion,:);
    smallFit = lsqnonneg(blobMat,profSmall(:));                            % fit to blobs using least-squares
     
    % embed into larger space
    smallFitBig           = zeros(se.movY*se.movX,1);
    smallFitBig(thesePix) = smallFit;
    blobFits(:,cc)        = smallFitBig;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

function tc = initTimeCourseCells(whichCells, nFrames)

% initialize storage of final results
tc.tcLSQ            = nan(nFrames,length(whichCells));                     % 
tc.tcSEUDO          = nan(nFrames,length(whichCells));                     % final output of the SEUDO-LS comparison
tc.lsqCosts         = nan(nFrames,length(whichCells));                     % Least-squares cost
tc.bobCosts         = nan(nFrames,length(whichCells));                     % Total cost of SEUDO cell construction
tc.tcCellsWithBlobs = nan(nFrames,length(whichCells));                     % Total cost for the known cells and the SEUDO cell construction
tc.tcBlobs          = cell(length(whichCells),1);                          %
tc.tcOtherCellsLSQ_all = cell(length(whichCells),1);                       % 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [md, whichProfiles, whichPixels, tcBlobs] = collectMetaInfo(se,...
                             whichCells, tcBlobs, one_blob, sigma2, params)

whichPixels      = false(se.movY*se.movX,length(whichCells));              % which pixels are used in each cell
whichProfiles    = false(se.nCells,length(whichCells));                    % which profiles to include when performing SEUDO on each cell
md               = cell(length(whichCells),1);                             % Initialize the cell array that holds all the meta data

for cc = 1:length(whichCells)
    
    thisCell  = whichCells(cc);
    prof      = se.profiles(:,:,thisCell);
    
    % identify pixel range for this cell
    % get outer bounds of area to which to add padding
    if params.useCOM                                                       % just the COM
        theCom = seudo.computeRoiCOMs(prof);
        ob = round(theCom([2 2 1 1]));
    else                                                                   % outer bound of non-zero pixels
        [~,ob] = seudo.computeRoiCOMs(prof);
    end
    
    % add padding
    ob = [max(1,ob(1)-params.padSpace) min(se.movY,ob(2)+params.padSpace)  max(1,ob(3)-params.padSpace) min(se.movX,ob(4)+params.padSpace)];
    
    % identify which subset of the original this is
    thesePix                          = false(se.movY,se.movX);
    thesePix(ob(1):ob(2),ob(3):ob(4)) = true;
    whichPixels(:,cc)                 = thesePix(:);
    
    if params.saveBlobTimeCourse
        tcBlobs{cc} = nan(se.movF,sum(whichPixels(:,cc)));            % initialize blob time courses if requested
    end

    % identify which cells to include along with this cell
    pixPerProfile = squeeze(sum(sum(se.profiles(ob(1):ob(2),ob(3):ob(4),:)>0)));
    whichProfiles(:,cc) = pixPerProfile > params.minPixForInclusion;
    
    % ensure this cell gets included
    whichProfiles(thisCell,cc) = true;
    
    md{cc}.cellIndexWithinTheseCells = find(find(whichProfiles(:,cc)) == thisCell);
    
    % Get Useful Sizes
    md{cc}.nY_all = length(ob(1):ob(2));                                   % Get size of this region
    md{cc}.nX_all = length(ob(3):ob(4));                                   % -----------------
    nY            = md{cc}.nY_all;
    nX            = md{cc}.nX_all;
    nBlobs        = nY*nX;  
    nCells        = sum(whichProfiles(:,cc));                                                          

    % reshape profiles into a matrix (each column one profile)
    md{cc}.rois_all = reshape(se.profiles(ob(1):ob(2),ob(3):ob(4),whichProfiles(:,cc)),nY*nX,nCells);
    rois            = md{cc}.rois_all;

    md{cc}.lambdas0_all = [params.lambdaProf*ones(1,nCells), ...
                                       params.lambdaBlob*ones(1,nBlobs)]'; % set lambda values
    
    % set constants automatically
    k1 = 2*sigma2*md{cc}.lambdas0_all;
    k2 = 2*sigma2*(log(1./params.p-1) - sum(log(md{cc}.lambdas0_all(md{cc}.lambdas0_all>0))));
    k3 = 2*sigma2*(log(1./params.p-1));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set solver options
    md{cc}.ls_opt_all            = struct;                                 % Create TFOCS options struct
    md{cc}.ls_opt_all.tol        = 0.01;                                   % Set algorithmic tolerance (for TFOCS)
    md{cc}.ls_opt_all.nonneg     = true;                                   % Set non-negative flag (for TFOCS)
    md{cc}.ls_opt_all.printEvery = 0;                                      % Suppress TFOCS output
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Create TFOCS Linear Operator
    
    md{cc}.prof_norms_all  = sqrt(sum(rois.^2,1));                         % Get only the profile norms
    md{cc}.normFactors_all = [md{cc}.prof_norms_all,ones(1,nBlobs)];       % Calculate the norms of all the factors used in the generative model
    md{cc}.lambdas_all     = k1./md{cc}.normFactors_all';                  % scale lambdas to match
    md{cc}.apply_blobs_all = @(z) reshape(convn(reshape(z, nY,nX),...
                                             one_blob, 'same'), nX*nY, 1); % Make an anonymous function that applies the blobs accross the image via convolution
    
    if strcmp(params.dyn_mod,'static')
%         Af   = @(z) rois*(vec(z(1:nCells))./(md{cc}.prof_norms_all.')) + reshape(conv_freq2(...
%          reshape(z(nCells+(1:nBlobs)),nY,nX),one_blob,'same'), nBlobs, 1); % Make forward operator
%         At   = @(z) [(rois.'*vec(z))./(md{cc}.prof_norms_all.'); reshape(conv_freq2(...
%                           reshape(z,nY,nX), one_blob, 'same'), nBlobs,1)]; % Make transpose operator
        Af   = @(z) rois*(vec(z(1:nCells))./(md{cc}.prof_norms_all.')) + reshape(convn(...
         reshape(z(nCells+(1:nBlobs)),nY,nX),one_blob,'same'), nBlobs, 1); % Make forward operator
        At   = @(z) [(rois.'*vec(z))./(md{cc}.prof_norms_all.'); reshape(convn(...
                          reshape(z,nY,nX), one_blob, 'same'), nBlobs,1)]; % Make transpose operator
        md{cc}.Xmat_all =  linop_handles([nBlobs, nCells+nBlobs], Af, At, 'R2R');           % Create a linear operator handle for TFOCS
    end
    md{cc}.ls_opt_all.L          = 2*linop_normest(md{cc}.Xmat_all).^2;

    switch params.dyn_mod
        case 'rwl1df'
            dp.xis    = dp.beta./md{cc}.lambdas_all(md{cc}.lambdas_all>0);
            dp.alphas = 0.5*dp.tau./dp.xis - 1;
            dp.k1     = k1;
            dp.k2     = k2;
            dp.k3     = k3;
        case 'bpdndf'
            dp.k1     = k1;
            dp.k2     = k2;
            dp.k3     = k3;
        otherwise
            dp.k1     = k1;
            dp.k2     = k2;
            dp.k3     = k3;
    end
    md{cc}.dp_all = dp;
    
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

function updateFrames = createFramesToUpdateOn(se, initialFrames, params)

% pick on which frames to show update
initialFrames = initialFrames(initialFrames<se.movF*.2);
prcFrames     = round(linspace(0,se.movF,21));
prcFrames     = prcFrames(prcFrames>initialFrames(end));
updateFrames  = [initialFrames prcFrames];
if ~params.verbose, updateFrames = []; end                                 %

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tcTemp = runOneSEUDO(isFirst, md, thisFrameFull, whichProfiles, whichPixels, prevTc, fitBOB, dyn_mod, optimizerOpt)

% How to call this function:
% tcTemp = runOneSEUDO(se, cc, ff==1, md{cc}, thisFrameFull, ...
%                          whichPixels(:,cc), whichPixels(:,cc), ...
%                               tc.tcSEUDO(ff-1,:).', tcBlobs{cc}(ff,:))

rois   = md.rois_all;                                                      % 
nY     = md.nY_all;                                                        % 
nX     = md.nX_all;                                                        % 
nCells = sum(whichProfiles);                                               % 

thisFrame       = thisFrameFull(whichPixels);                              % Extract the relevant portion of the current frame
tcLSQ_thisFrame = ((rois'*rois)\(rois'*thisFrame))';                       % least squares solution
fitWeights      = md.Xmat_all(thisFrame,2)-md.lambdas_all;                 % Initialize RWL1-DF estimate to zero
fitWeights(fitWeights<0) = 0;
if (~isFirst)&&(strcmp(dyn_mod,'bpdndf')||strcmp(dyn_mod,'rwl1df'))
    tc_past = [prevTc;vec(fitBOB)];
    md.dp_all.k2 = update_addconst(vec(fitBOB),params,md.lambdas0_all);
else; tc_past = fitWeights;
end
switch dyn_mod
    case 'rwl1df'
        fitWeights = singleRWL1DFStep(thisFrame, tc_past, ...
            fitWeights, one_blob, rois, md.lambdas_all, ...
            md.prof_norms_all, [nX,nY,nCells], md.ls_opt_all, md.dp_all); % Update weights using RWL1DF
    case 'bpdndf'
        fitWeights = singleBPDNDFStep(thisFrame, tc_past, ...
            fitWeights, one_blob, rois, md.lambdas_all, md.dp_all, ...
            md.prof_norms_all, [nX,nY,nCells], md.ls_opt_all);
    otherwise
        if strcmp(optimizerOpt,'TFOCS')
            md.ls_opt_all = rmfield(md.ls_opt_all,'L');
            fitWeights = solver_L1RLS(md.Xmat_all, thisFrame, md.lambdas_all, fitWeights , md.ls_opt_all);    % get best fit weights
        elseif strcmp(optimizerOpt,'prox')
            fitWeights = proxGradDescentNonnegWL1(md.Xmat_all, thisFrame, tc_past, md.lambdas_all, md.ls_opt_all);
        end
end

tcTemp.fitWeights = fitWeights;
fitWeights = fitWeights./md.normFactors_all';                              % undo the norm factors that were needed to condition A
fitX       = fitWeights(1:nCells);                                         % separate weight vector into found cell activation values (X) ...
fitBOB     = fitWeights(nCells+1:end);                                     % ...and not-found cell activation values (Bunch Of Blobs)

% decide whether this frame is best fit by LSQ or the cells plus blobs
lsqCost = norm(rois*tcLSQ_thisFrame'-thisFrame,2)^2;                       % compute the least-squares cost of each fit
if (~isFirst)&&(strcmp(dyn_mod,'bpdndf')||strcmp(dyn_mod,'rwl1df'))
    bobCost = dynamicSEUDOcost(fitWeights, thisFrame, tc_past, rois, ...
                           md.lambdas_all, md.dp_all, md.apply_blobs_all); % calculate the dynamic SEUDO cost
else
    bobCost = norm(thisFrame-rois*fitX-md.apply_blobs_all(fitBOB),2)^2 ...
                                + sum(abs(md.dp_all.k1(:).*fitWeights(:)))-md.dp_all.k2; % Cost w/ blobs is a log-likelihood - in this case the LASSO cost [ k2 = 2*sigma2*(log(1./q-1) + sum(log(lambdas0((lambdas0>0))))) ]
end

% if lsqCost<bobCost; fprintf('!!!!!!!!!\n'); end %% REMOVE ME

% choose the fit with lower cost
if lsqCost<bobCost*(-Inf)                                                  % IF LEAST-SQUARES HAS A LOWER COST
    fitFancy = tcLSQ_thisFrame;                                            % Save the least-squares solution
    fitBOB   = 0*fitBOB;                                                   % Set the blob costs back to zero
    fitX     = 0*fitX;                                                     % Set the blob costs back to zero
else                                                                       % IF THE SEUDO MODEL HAS A LOWER COST
    fitFancy = fitX;                                                       % Save the SEUDO model solution
end

% store fits in time course variables
tcTemp.tcSEUDO             = fitFancy(md.cellIndexWithinTheseCells);       %
tcTemp.tcCellsWithBlobs    = fitX(md.cellIndexWithinTheseCells);           %
tcTemp.tcLSQ               = tcLSQ_thisFrame(md.cellIndexWithinTheseCells);%
tcTemp.tcLSQ_thisFrame     = tcLSQ_thisFrame;                              %
tcTemp.fitBOB              = fitBOB;                                       % save blobs
tcTemp.lsqCosts            = lsqCost;                                      % store the Least-squares costs
tcTemp.bobCosts            = bobCost;                                      % store the blob costs
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function bobCost = dynamicSEUDOcost(fitWeights, dF_now, tc_old, rois, lambdas, dp, apply_blobs)

% bobCost = dynamicSEUDOcost(fitWeights, dF_now, tc_old, rois, lambdas, dp, p, apply_blobs)
%
%
%
% 2018 - Adam Charles

% Some initializations

nRois = numel(fitWeights)-numel(dF_now);

fitX   = fitWeights(1:nRois);
fitBOB = fitWeights(nRois+1:end);

% Calculate cost

switch dyn_mod
    case 'rwl1df'
        alphas  = dp.alphas;
        eta     = dp.eta;
        xis     = dp.xis;
        k1      = dp.k1;
        k3      = dp.k3;
        lambdas = lambdas(nRois+1:end);
        bobCost = norm(dF_now-rois*fitX-apply_blobs(fitBOB),2)^2 -k3 ...
           + sum(k1(nRois+1:end).*(alphas(nRois+1:end).*log((abs(tc_old(nRois+1:end)) + eta)./xis(nRois+1:end))...
           + alphas(nRois+1:end).*log(lambdas(lambdas>0).*fitWeights((nRois+1):end) ...
               + (abs(tc_old(nRois+1:end)) + eta)./xis(nRois+1:end))...
                                  + log(2) + alphas(nRois+1:end).*log(alphas(nRois+1:end))));% Cost w/ blobs is a log-likelihood - in this case the LASSO cost [ k2 = 2*sigma2*(log(1./q-1) + sum(log(lambdas0((lambdas0>0))))) ]
    case 'bpdndf'
        k1      = dp.k1;
        k2      = dp.k2;
        kappa   = dp.kappa;
        bobCost = norm(dF_now-rois*fitX-apply_blobs(fitBOB),2)^2 ...
                        + kappa*norm(vec(fitWeights) - vec(tc_old))^2 ...
                                    + sum(abs(k1(:).*fitWeights(:))) - k2; % Cost w/ blobs and dynamics is ...
    otherwise
        error('dynamics type %s not recoanized',dyn_mod)
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function k2 = update_addconst(fitBoB, params, lambdas0)

% k2 = update_addconst(fitBoB, params, lambdas0)
% 
% Code to update the additive constant for the maximum likelihood
% comparison in the SEUDO code. This is especially necessary in the dynamic
% SEUDO code when the constants change at each iteration. 
% 
% 2018 - Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the offset constant

switch params.dyn_mod
    case 'rwl1df'
        k2 = 2*params.sigma2*(log(1./params.p-1) + sum(log(1./lambdas0((lambdas0>0)))));
    case 'bpdndf'  
        logZ = params.lambdaBlob*sum(abs(fitBoB)) ...
            - params.lambdaBlob*numel(fitBoB)./(4*params.kappa) ...
            - sum(log(normcdf(sqrt(2*params.kappa)*fitBoB - ...
                    params.lambdaBlob/sqrt(2*params.kappa))));
        k2 = 2 *params.sigma2*(log(1./params.p-1) - logZ );
    otherwise
        k2 = 2*params.sigma2*(log(1./params.p-1) - sum(log(lambdas0(lambdas0>0))));
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fitWeights = singleRWL1DFStep(dF_now, tc_old, fitWeights, one_blob, rois, lambdas, prof_norms, szs, ls_opt, dp)

% fitWeights = singleRWL1DFStep(dF_now, tc_old, fitWeights, one_blob, rois, lambdas, prof_norms, szs, ls_opt)
%
% Function that calculates a single step using the re-weighted l1 dynamic
% filtering regularization. 
%
% 2018 - Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nX     = szs(1);
nY     = szs(2);
nRois  = szs(3);
nBlobs = nX*nY;

eta  = dp.eta;
tau  = dp.tau;
beta = dp.beta;

if (dp.n_rw_iter <= 0)&&(~isempty(tc_old))
    lambdas_rw = tau./(eta+abs(tc_old));                                   % Initialize the lambda values (weighted LASSO parameters) for the RWL1 iterations
else
    lambdas_rw = ones(size(lambdas));                                      % Initialize the lambda values (weighted LASSO parameters) for the RWL1 iterations
end
if ~isempty(tc_old)
    for ll = 1:dp.n_rw_iter
        Af   = @(z) (rois*diag(1./prof_norms))*(vec(z(1:nRois))./lambdas_rw(1:nRois)) + ...
                 reshape(convn(reshape(z(nRois+(1:nX*nY))./lambdas_rw(nRois+(1:nX*nY)), ...
                                     nY,nX), one_blob, 'same'), nX*nY, 1); % Make forward operator
        At   = @(z) [(rois*diag(1./prof_norms)).'*vec(z); ...
                        reshape(convn(reshape(z, nY,nX), one_blob, ...
                                          'same'), nBlobs,1)]./lambdas_rw; % Make transpose operator
        Xmat_all =  linop_handles([nX*nY, nRois+nBlobs], Af, At, 'R2R');   % Create a linear operator handle for TFOCS

        fitWeights = solver_L1RLS(Xmat_all,dF_now,lambdas,...
                                                      fitWeights,ls_opt);  % get best fit estimate
        lambdas_rw = tau./(eta+beta*abs(tc_old) + ...
                                               (1-beta)*abs(fitWeights));  % Update the weighted LASSO parameters based on both the current estimate and the previous estimate
    end
else
    Af   = @(z) (rois*diag(1./prof_norms))*(vec(z(1:nRois))./lambdas_rw(1:nRois)) + ...
             reshape(convn(reshape(z(nRois+(1:nX*nY))./lambdas_rw(nRois+(1:nX*nY)), ...
                                     nY,nX), one_blob, 'same'), nX*nY, 1); % Make forward operator
    At   = @(z) [(rois*diag(1./prof_norms)).'*vec(z); ...
                    reshape(convn(reshape(z, nY,nX), one_blob, ...
                                          'same'), nBlobs,1)]./lambdas_rw; % Make transpose operator
    Xmat_all  = linop_handles([nX*nY, nRois+nBlobs], Af, At, 'R2R');       % Create a linear operator handle for TFOCS
    fitWeights = solver_L1RLS(Xmat_all,dF_now,lambdas,...
                                                      fitWeights,ls_opt);  % get best fit estimate
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fitWeights = singleBPDNDFStep(dF_now, tc_old, fitWeights, one_blob, rois, lambdas, dp, prof_norms, szs, ls_opt)

%
%
%
%
% 2018 - Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up sized and extract parameters

nX        = szs(1);
nY        = szs(2);
nRois     = szs(3);
nBlobs    = nX*nY;
kappa     = dp.kappa;
kappaProf = dp.kappaProf;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the activity estimate using BPDN-DF

if ~isempty(tc_old)                                                        % Run one BPDN-DF iteration
    c1 = sqrt(2*kappaProf*dp.sigma2);
    c2 = sqrt(2*kappa*dp.sigma2);
    k  = [c1*ones(nRois,1);c2*ones(nBlobs,1)];                             % Set up the necessary constants
    y  = [dF_now; k.*tc_old];                                              % Augment the least-squares "measurement" vector with the history portion
    
    Af   = @(z) [(rois*diag(1./prof_norms))*vec(z(1:nRois)) + ...
                 reshape(convn(reshape(z(nRois+(1:nX*nY)),nY,nX), ...
                                       one_blob, 'same'), nX*nY, 1); k.*z]; % Make forward operator
    At   = @(z) [(rois*diag(1./prof_norms)).'*vec(z(1:nX*nY)); ...
              reshape(convn(reshape(z(1:nX*nY), nY,nX), ...
                       one_blob, 'same'), nBlobs,1)] + k.*z((nX*nY+1):end); % Make transpose operator
    Xmat =  linop_handles([nX*nY+nRois+nBlobs, nRois+nBlobs],Af,At,'R2R'); % Create a linear operator handle for TFOCS

    fitWeights = solver_L1RLS(Xmat,y,lambdas,fitWeights,ls_opt);           % get best fit estimate
else                                                                       % On the first iteration there is no history to use
    Af   = @(z) (rois*diag(1./prof_norms))*vec(z(1:nRois)) + reshape(...
             convn(reshape(z(nRois+(1:nX*nY)), nY,nX), one_blob, ...
                                                       'same'), nX*nY, 1); % Make forward operator
    At   = @(z) [(rois*diag(1./prof_norms)).'*vec(z); reshape(convn(...
                         reshape(z, nY,nX), one_blob, 'same'), nBlobs,1)]; % Make transpose operator
    Xmat  = linop_handles([nX*nY, nRois+nBlobs], Af, At, 'R2R');           % Create a linear operator handle for TFOCS
    fitWeights = solver_L1RLS(Xmat,dF_now,lambdas,fitWeights,ls_opt);      % get best fit estimate
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% END OF FILE
