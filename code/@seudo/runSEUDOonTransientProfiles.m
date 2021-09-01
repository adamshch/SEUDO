function se = runSEUDOonTransientProfiles(se,varargin)
% seudoObject.runSEUDOonTransientProfiles(...)
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
% 2019 Adam Charles

if ~exist('tfocs','file')
    error('TFOCS must be installed to use SEUDO, download from http://cvxr.com/tfocs/download/')
end
CODE_VERSION = 1;                                                          % Set code version

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input Parameter parsing

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

params.minPixForInclusion = max(params.minPixForInclusion,1);              % ensure at least 1
tc       = initTimeCourseCells(se, whichCells);                            % initialize storage of final results
tcBlobs  = cell(length(whichCells),1);                                     % initialize the storage of the blob timecourses
one_blob = makeOneBlob(params.blobRadius);                                 % make one blob

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% initialize storage of cell-specific variables (1 cell index = 1 cell)

[md, whichProfiles, whichPixels, tc.tcBlobs] = collectMetaInfo(se, ...
                 whichCells, tc.tcBlobs, one_blob, params.sigma2, params); % Get a bunch of meta data containing the locations/pixels/parameters for SEUDO for each cell 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get blob fit to each cell

blobFits = fitBlobsToProfiles(se, whichCells, params.padSpace, params.blobRadius);  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run per-transient analysis

startTime = now;
fprintf('Running')
for cc = 1:length(whichCells)                                              % Iterate and run SEUDO on each cell
    % get cell-specific info
    rois   = md.rois_all{cc};                                              % 
    nY     = md.nY_all{cc};                                                % 
    nX     = md.nX_all{cc};                                                % 
    nCells = sum(whichProfiles(:,cc));                                     % 
    dp     = md.dp_all{cc};                                                %     
    
    for ff=1:md.nTran(cc)       
        % Iterate over frames
        thisFrame = getCurrentData(se, 'transProfile', cc, ff);            % Pick out the current frame
        
        tcLSQ_thisFrame = ((rois'*rois)\(rois'*thisFrame))';               % least squares solution for this frame   
        
        %% perform SEUDO on this frame
        fitWeights = zeros(nCells+nY*nX,1);                                    % Initialize RWL1-DF estimate to zero
        fitWeights = solver_L1RLS(md.Xmat_all{cc}, thisFrame, ...
                      md.lambdas_all{cc}, fitWeights , md.ls_opt_all{cc}); % get best fit weights
        % NOTE: In order to NOT use TFOCS, uncomment this line and block above to make H_all etc.
%                 fitWeights = quadprog(H_all{cc}, (h_all{cc}*thisFrame+lambdas_all{cc}), [], [], [], [], zeros(nX*nY+nCells,1), Inf*ones(nX*nY+nCells,1), fitWeights, optimset('Display', 'off','TolFun',1e-14,'Algorithm','trust-region-reflective')); 

        fitWeights = fitWeights./md.normFactors_all{cc}';                  % undo the norm factors that were needed to condition A
        fitX       = fitWeights(1:nCells);                                 % separate weight vector into found cell activation values (X) and
        fitBOB     = fitWeights(nCells+1:end);                             % ... and the not-found cell activation values (Bunch Of Blobs)
        
        %% decide whether this frame is best fit by LSQ or the cells plus blobs
        
        lsqCost = norm(rois*tcLSQ_thisFrame'-thisFrame,2)^2;               % compute the least-squares cost of each fit
        bobCost = norm(thisFrame-rois*fitX-md.apply_blobs_all{cc}(fitBOB),2)^2 ...
                                + sum(abs(dp.k1(:).*fitWeights(:)))-dp.k2; % Cost w/ blobs is a log-likelihood - in this case the LASSO cost [ k2 = 2*sigma2*(log(1./q-1) + sum(log(lambdas0((lambdas0>0))))) ]

        % choose the fit with lower cost
        if lsqCost<bobCost                                                 % IF LEAST-SQUARES HAS A LOWER COST
            fitFancy = tcLSQ_thisFrame;                                    % Save the least-squares solution
            fitBOB   = 0*fitBOB;                                           % Set the blob costs back to zero
            fitX     = 0*fitX;                                             % Set the blob costs back to zero
        else                                                               % IF THE SEUDO MODEL HAS A LOWER COST
            fitFancy = fitX;                                               % Save the SEUDO model solution
        end
        
        % store fits in time course variables
        tc.tcSEUDO{cc}(ff)          = fitFancy(md.cellIndexWithinTheseCells{cc});
        tc.tcCellsWithBlobs{cc}(ff) = fitX(md.cellIndexWithinTheseCells{cc});
        tc.tcLSQ{cc}(ff)            = tcLSQ_thisFrame(md.cellIndexWithinTheseCells{cc});
        
        if params.saveOtherCellTimeCourses
            tc.tcOtherCellsLSQ_all{cc}(ff,:) = tcLSQ_thisFrame;
        end
        
        % save blobs
        if params.saveBlobTimeCourse
            tc.tcBlobs{cc}(ff,:) = fitBOB;
        else
            tcBlobs{cc}(ff,1)  = sum(fitBOB);
            prof               = se.profiles(:,:,cc) > 0;
            tcBlobs{cc}(ff,2)  = sum(fitBOB(prof(whichPixels(:,cc))));
        end
        
        tc.lsqCosts{cc}(ff) = lsqCost;                                     % store the Least-squares costs
        tc.bobCosts{cc}(ff) = bobCost;                                     % store the blob costs
        md.dp_all{cc}       = dp;                                          % store updated dp
    end
    fprintf('.')
end
fprintf('\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output parsing

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
se.tcSeudo.perTran.tc     = tc.tcSEUDO;                                    % SEUDO time courses 
se.tcSeudo.perTran.tcLSQ  = tc.tcLSQ;                                      % LSQ time courses 
se.tcSeudo.perTran.params = params;                                        % parameters
se.tcSeudo.perTran.extras = extras;                                        % intermediate computations for debugging

if params.verbose
    fprintf('   completed in %0.1f min (%s), results stored in field perTran\n',...
        (now-startTime)*24*60,datestr(now))
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

p.addParameter('whichCells'               , []);
p.addParameter('saveBlobTimeCourse'       , false);
p.addParameter('saveOtherCellTimeCourses' , true);

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

function tc = initTimeCourseCells(se, whichCells)

% initialize storage of final results
tc.tcLSQ               = cell(length(whichCells),1);                       %
tc.tcSEUDO             = cell(length(whichCells),1);                       % final output of the SEUDO-LS comparison
tc.lsqCosts            = cell(length(whichCells),1);                       % Least-squares cost
tc.bobCosts            = cell(length(whichCells),1);                       % Total cost of SEUDO cell construction
tc.tcCellsWithBlobs    = cell(length(whichCells),1);                       % Total cost for the known cells and the SEUDO cell construction
tc.tcBlobs             = cell(length(whichCells),1);                       %
tc.tcOtherCellsLSQ_all = cell(length(whichCells),1);                       % 

for cc = 1:length(whichCells)
    thisCell  = whichCells(cc);
    nTran     = size(se.tcDefault.transientInfo(thisCell).shapes,3);
    tc.tcLSQ{cc}               = nan(nTran,1);                             %
    tc.tcSEUDO{cc}             = nan(nTran,1);                             % final output of the SEUDO-LS comparison
    tc.lsqCosts{cc}            = nan(nTran,1);                             % Least-squares cost
    tc.bobCosts{cc}            = nan(nTran,1);                             % Total cost of SEUDO cell construction
    tc.tcCellsWithBlobs{cc}    = nan(nTran,1);                             % Total cost for the known cells and the SEUDO cell construction
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [md, whichProfiles, whichPixels, tcBlobs] = collectMetaInfo(se,...
                             whichCells, tcBlobs, one_blob, sigma2, params)

whichPixels      = false(se.movY*se.movX,length(whichCells));              % which pixels are used in each cell
whichProfiles    = false(se.nCells,length(whichCells));                    % which profiles to include when performing SEUDO on each cell
 
md.nX_all          = cell(length(whichCells),1);                           % 
md.nY_all          = cell(length(whichCells),1);                           % 
md.rois_all        = cell(length(whichCells),1);                           % 
md.prof_norms_all  = cell(length(whichCells),1);                           % 
md.normFactors_all = cell(length(whichCells),1);                           % 
md.lambdas_all     = cell(length(whichCells),1);                           % vector of lambda values for all components 
md.apply_blobs_all = cell(length(whichCells),1);                           % anonymous function that convolves an image with the Gaussian kernel
md.dp_all          = cell(length(whichCells),1);                           % 
md.lambdas0_all    = cell(length(whichCells),1);                           % 
md.ls_opt_all      = cell(length(whichCells),1);                           % 
md.Xmat_all        = cell(length(whichCells),1);                           % 
md.nTran           = zeros(length(whichCells),1);                          % Vector of number of transients per cell

md.cellIndexWithinTheseCells = cell(length(whichCells),1);                 % 

for cc = 1:length(whichCells)
    
    thisCell  = whichCells(cc);
    prof      = se.profiles(:,:,thisCell);
    md.nTran(cc) = size(se.tcDefault.transientInfo(thisCell).shapes,3);
    
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
        tcBlobs{cc} = nan(md.nTran(cc),sum(whichPixels(:,cc)));            % initialize blob time courses if requested
    end

    % identify which cells to include along with this cell
    pixPerProfile = squeeze(sum(sum(se.profiles(ob(1):ob(2),ob(3):ob(4),:)>0)));
    whichProfiles(:,cc) = pixPerProfile > params.minPixForInclusion;
    
    % ensure this cell gets included
    whichProfiles(thisCell,cc) = true;
    
    md.cellIndexWithinTheseCells{cc} = find(find(whichProfiles(:,cc)) == thisCell);
    
    % Get Useful Sizes
    md.nY_all{cc} = length(ob(1):ob(2));                                   % Get size of this region
    md.nX_all{cc} = length(ob(3):ob(4));                                   % -----------------
    nY = md.nY_all{cc};
    nX = md.nX_all{cc};
    nBlobs = nY*nX;  
    nCells = sum(whichProfiles(:,cc));                                                          

    % reshape profiles into a matrix (each column one profile)
    md.rois_all{cc} = reshape(se.profiles(ob(1):ob(2),ob(3):ob(4),whichProfiles(:,cc)),nY*nX,nCells);
    rois            = md.rois_all{cc};

    md.lambdas0_all{cc} = [params.lambdaProf*ones(1,nCells), ...
                                       params.lambdaBlob*ones(1,nBlobs)]'; % set lambda values
    
    % set constants automatically
    k1 = 2*sigma2*md.lambdas0_all{cc};
    k2 = 2*sigma2*(log(1./params.p-1) - sum(log(md.lambdas0_all{cc}(md.lambdas0_all{cc}>0))));
    k3 = 2*sigma2*(log(1./params.p-1));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set solver options
    md.ls_opt_all{cc}            = struct;                                                % Create TFOCS options struct
    md.ls_opt_all{cc}.tol        = 0.01;                                                  % Set algorithmic tolerance (for TFOCS)
    md.ls_opt_all{cc}.nonneg     = true;                                                  % Set non-negative flag (for TFOCS)
    md.ls_opt_all{cc}.printEvery = 0;                                                     % Suppress TFOCS output
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Create TFOCS Linear Operator
    
    md.prof_norms_all{cc}  = sqrt(sum(rois.^2,1));                            % Get only the profile norms
    md.normFactors_all{cc} = [md.prof_norms_all{cc},ones(1,nBlobs)];             % Calculate the norms of all the factors used in the generative model
    md.lambdas_all{cc}     = k1./md.normFactors_all{cc}';                        % scale lambdas to match
    md.apply_blobs_all{cc} = @(z) reshape(convn(reshape(z, nY,nX),...
                                             one_blob, 'same'), nX*nY, 1); % Make an anonymous function that applies the blobs accross the image via convolution
    
    if strcmp(params.dyn_mod,'static')
        Af   = @(z) (rois*diag(1./md.prof_norms_all{cc}))*vec(z(1:nCells)) + reshape(convn(...
            reshape(z(nCells+(1:nX*nY)), nY,nX), one_blob, 'same'), nX*nY, 1); % Make forward operator
        At   = @(z) [(rois*diag(1./md.prof_norms_all{cc})).'*vec(z); reshape(convn(...
            reshape(z, nY,nX), one_blob, 'same'), nBlobs,1)]; % Make transpose operator
        md.Xmat_all{cc} =  linop_handles([nX*nY, nCells+nBlobs], Af, At, 'R2R');           % Create a linear operator handle for TFOCS
    end
    
    % NOTE: In order to NOT use TFOCS, uncomment these lines and the switch below (b/w TFOCS & quadprog)
%     miniBlobsMat_all{cc} = makeBlobsMatrix(nY,nX,blurRad,clipHeight);      % Create the matrix of blobs
%     H_all{cc} = 2*([rois*diag(1./prof_norms_all{cc}), miniBlobsMat_all{cc}].'*[rois*diag(1./prof_norms_all{cc}), miniBlobsMat_all{cc}]);
%     H_all{cc} = sparse(0.5*(H_all{cc}+H_all{cc}'));
%     h_all{cc} = -2*[rois*diag(1./prof_norms_all{cc}), miniBlobsMat_all{cc}].';

    dp.k1      = k1;
    dp.k2      = k2;
    dp.k3      = k3;
    md.dp_all{cc} = dp;
    
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [thisFrame,varargout] = getCurrentData(se, dataType, cc, ff, varargin)

if nargin > 4
    whichPixels = varargin{1};
else
    whichPixels = [];
end

if nargin > 5
    frameSlidingWindow = varargin{1};
else
    frameSlidingWindow = [];
end

if strcmp(dataType,'transProfile')
    thisFrame  = reshape(se.tcDefault.transientInfo(cc).shapes(:,:,ff),[],1);% Get current frame
elseif strcmp(dataType,'frame')
    currFrame = reshape(se.getFrame(ff),[],1);                             % get current frame
    if ff <= dsT
        frameSlidingWindow(:,ff) = currFrame;                              % the first few times, build up a matrix of the initial frames
    else
        frameSlidingWindow = [frameSlidingWindow(:,2:end) currFrame];      % subequent times, delete the oldest frame and appent the current frame
    end
    thisFrameFull = nanmean(frameSlidingWindow,2);                         % average the last few frames

    thisFrame       = thisFrameFull(whichPixels(:,cc));                    % Get current frame
else
end

if nargout > 1
    varargout{1} = frameSlidingWindow;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






