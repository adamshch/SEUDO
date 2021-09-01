function se = estimateTimeCoursesWithSEUDO_perCellParams(se,varargin)
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
% 2018 Jeff Gauthier



CODE_VERSION = 1;

% Get Input Parameters

% parase parameters
p = inputParser;

p.addParameter('p'                  , 0.00001, @(x)isnumeric(x) & all( x(:)>=0 & x(:)<=1 ));   
p.addParameter('sigma2'             , .01, @isnumeric);   
p.addParameter('lambdaBlob'         , 20, @isnumeric); 
p.addParameter('blobRadius'         , 1.2, @isnumeric);

p.addParameter('dsTime'             , 1, @isnumeric);
p.addParameter('lambdaProf'         , 0, @isnumeric);
p.addParameter('minPixForInclusion' , 1, @isnumeric);
p.addParameter('padSpace'           , 10, @isnumeric);
p.addParameter('useCOM'             , false, @islogical);

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
p.addParameter('saveBlobTimeCourse'       , false, @islogical);
p.addParameter('saveOtherCellTimeCourses' , true, @islogical);



% just validate inputs, if desired
q = inputParser;
q.KeepUnmatched = true;
q.addParameter('validateParams' , false, @(x)isscalar(x)&&(islogical(x) || ismember(x,[0 1])));
parse(q,varargin{:});
if q.Results.validateParams
    % if inputs can be parsed successfully, return empty
    % otherwise return the error
    se = [];
    try parse(p,varargin{:}); catch err, se = err; end
    return
end


% just show defaults, if desired
q = inputParser;
q.KeepUnmatched = true;
q.addParameter('showDefaults' , false, @(x)isscalar(x)&&(islogical(x) || ismember(x,[0 1])));
parse(q,varargin{:});
if q.Results.showDefaults
    % get default values
    parse(p);
    params = p.Results;
    params.version = CODE_VERSION;
    % if there are output variables, return struct
    if nargout > 0
        se = params;
        return
    else
        % otherwise, just display default values
        fprintf('\ndefault values for SEUDO:\n\n')
        disp(params)
        clear se
        return
    end
end


% parse inputs
parse(p,varargin{:});




% copy to struct that will be output
params = p.Results;

% set code version
params.version = CODE_VERSION;




% split into per-cell parameters

% get list
perCellParams = {'p','sigma2','lambdaBlob','blobRadius','lambdaProf','minPixForInclusion','padSpace','useCom'};

% verify each provided parameter has the correct number of elements
for pp = 1:length(perCellParams)
    if ~ismember(length(params.(perCellParams{pp})),[1 se.nCells])
        error('parameter ''%s'' must be a scalar or %d-length vector',...
            perCellParams{pp},se.nCells)
    end
end




% ESTIMATE TIME COURSES 

if isempty(params.whichCells)
    whichCells = 1:size(se.profiles,3);
else
    whichCells = params.whichCells;
end



% put parameters into convenient names
padSpace = params.padSpace;



% sizes to keep track of
nFrames = se.movF;                                                         % Number of time-frames
nCellsTotal = length(whichCells);                                          % Number of cells


% initialize storage of final results
tcLSQ            = nan(nFrames,nCellsTotal);
tcSEUDO          = nan(nFrames,nCellsTotal);                               % final output of the SEUDO-LS comparison
lsqCosts         = nan(nFrames,nCellsTotal);                               % Least-squares cost
bobCosts         = nan(nFrames,nCellsTotal);                               % Total cost of SEUDO cell construction
tcCellsWithBlobs = nan(nFrames,nCellsTotal);                               % Total cost for the known cells and the SEUDO cell construction
whichPixels      = false(se.movY*se.movX,nCellsTotal);                     % which pixels are used in each cell
whichProfiles    = false(se.nCells,length(whichCells));                    % which profiles to include when performing SEUDO on each cell


if params.saveBlobTimeCourse
   %tcBlobs      = nan(nFrames,se.movY*se.movX,nCellsTotal);                     % save time course of each blob
    tcBlobs      = cell(nCellsTotal,1);
else
    tcBlobs      = nan(nFrames,2,nCellsTotal);                             % save summed blob time course
    % tcBlobs(:,1,cellID) is all blobs
    % tcBlobs(:,2,cellID) is only blobs overlapping the profile
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set some parameters

dsT    = params.dsTime;                                                    % temporal downsample
p      = params.p;                                                         % p = prob(not found cell is active)
sigma2 = params.sigma2/dsT;                                                % noise variance
params.minPixForInclusion = max(params.minPixForInclusion,1);              % ensure at least 1


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make one blob

% parameters
blurRad     = params.blobRadius;
clipHeight  = 0.01;
blobsMatrix = makeBlobsMatrix(se.movY,se.movX,blurRad,clipHeight);         % Create the matrix of blobs
one_blob    = reshape(blobsMatrix(:,sub2ind([se.movY se.movX], ...            
               round(0.5*se.movY),round(0.5*se.movX))),se.movY,se.movX).'; % Isolate a single blob - start with a single frame with the blob in the middle ...
[imax,jmax] = find(one_blob==max(one_blob(:)),1);                          % ... Find the blob center ...
cropRad     = ceil(blurRad*2.5+eps);
one_blob    = one_blob(intersect(imax+(-cropRad:cropRad),1:se.movY),...       
                            intersect(jmax+(-cropRad:cropRad),1:se.movX)); % ... And crop the blob
one_blob    = one_blob/sqrt(sum(one_blob(:).^2));                          % Normalize the blob




% initialize storage of cell-specific variables (1 cell index = 1 cell)
nX_all          = cell(length(whichCells),1);                              % 
nY_all          = cell(length(whichCells),1);                              % 
rois_all        = cell(length(whichCells),1);                              % 
prof_norms_all  = cell(length(whichCells),1);                              % 
normFactors_all = cell(length(whichCells),1);                              % 
lambdas_all     = cell(length(whichCells),1);                              % vector of lambda values for all components 
apply_blobs_all = cell(length(whichCells),1);                              % anonymous function that convolves an image with the Gaussian kernel
dp_all          = cell(length(whichCells),1);                              % 
lambdas0_all    = cell(length(whichCells),1);                              % 
ls_opt_all      = cell(length(whichCells),1);                              % 
Xmat_all        = cell(length(whichCells),1);                              % 

cellIndexWithinTheseCells = cell(length(whichCells),1);                    % 
tcOtherCellsLSQ_all       = cell(length(whichCells),1);                    % 
    


for cc = 1:length(whichCells)
    
    thisCell = whichCells(cc);
    prof = se.profiles(:,:,thisCell);
    
    
    % identify pixel range for this cell
    
    % get outer bounds of area to which to add padding
    if params.useCOM % just the COM
        theCom = seudo.computeRoiCOMs(prof);
        ob = round(theCom([2 2 1 1]));
    else % outer bound of non-zero pixels
        [~,ob] = seudo.computeRoiCOMs(prof);
    end
    
    % add padding
    ob = [max(1,ob(1)-padSpace) min(se.movY,ob(2)+padSpace)  max(1,ob(3)-padSpace) min(se.movX,ob(4)+padSpace)];
    
    % identify which subset of the original this is
    thesePix = false(se.movY,se.movX);
    thesePix(ob(1):ob(2),ob(3):ob(4)) = true;
    whichPixels(:,cc) = thesePix(:);
    
    
    % initialize blob time courses
    if params.saveBlobTimeCourse
        tcBlobs{cc} = nan(nFrames,sum(whichPixels(:,cc)));
    end
    
    
    % identify which cells to include along with this cell
    pixPerProfile = squeeze(sum(sum(se.profiles(ob(1):ob(2),ob(3):ob(4),:)>0)));
    whichProfiles(:,cc) = pixPerProfile > params.minPixForInclusion;
    
    % ensure this cell gets included
    whichProfiles(thisCell,cc) = true;
    
    cellIndexWithinTheseCells{cc} = find(find(whichProfiles(:,cc)) == thisCell);
    
    
    % Get Useful Sizes
    nY_all{cc} = length(ob(1):ob(2));                                                  % Get size of this region
    nX_all{cc} = length(ob(3):ob(4));                                                  % -----------------
    nY = nY_all{cc};
    nX = nX_all{cc};
    nBlobs = nY*nX;  
    nCells = sum(whichProfiles(:,cc));                                                          
    

    
    % reshape profiles into a matrix (each column one profile)
    rois_all{cc} = reshape(se.profiles(ob(1):ob(2),ob(3):ob(4),whichProfiles(:,cc)),nY*nX,nCells);
    
    rois = rois_all{cc};
    
    
    lambdas0_all{cc} = [params.lambdaProf*ones(1,nCells), ...
        params.lambdaBlob*ones(1,nBlobs)]';  % set lambda values
    
    % set constants automatically
    k1 = 2*sigma2*lambdas0_all{cc};
    k2 = 2*sigma2*(log(1./p-1) - sum(log(lambdas0_all{cc}(lambdas0_all{cc}>0))));
    k3 = 2*sigma2*(log(1./p-1));
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set solver options
    ls_opt_all{cc}            = struct;                                                % Create TFOCS options struct
    ls_opt_all{cc}.tol        = 0.01;                                                  % Set algorithmic tolerance (for TFOCS)
    ls_opt_all{cc}.nonneg     = true;                                                  % Set non-negative flag (for TFOCS)
    ls_opt_all{cc}.printEvery = 0;                                                     % Suppress TFOCS output
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Create TFOCS Linear Operator
    
    prof_norms_all{cc}  = sqrt(sum(rois.^2,1));                            % Get only the profile norms
    normFactors_all{cc} = [prof_norms_all{cc},ones(1,nBlobs)];             % Calculate the norms of all the factors used in the generative model
    lambdas_all{cc}     = k1./normFactors_all{cc}';                        % scale lambdas to match
    apply_blobs_all{cc} = @(z) reshape(convn(reshape(z, nY,nX),...
                                             one_blob, 'same'), nX*nY, 1); % Make an anonymous function that applies the blobs accross the image via convolution
    
    if strcmp(params.dyn_mod,'static')
        Af   = @(z) (rois*diag(1./prof_norms_all{cc}))*vec(z(1:nCells)) + reshape(convn(...
            reshape(z(nCells+(1:nX*nY)), nY,nX), one_blob, 'same'), nX*nY, 1); % Make forward operator
        At   = @(z) [(rois*diag(1./prof_norms_all{cc})).'*vec(z); reshape(convn(...
            reshape(z, nY,nX), one_blob, 'same'), nBlobs,1)]; % Make transpose operator
        Xmat_all{cc} =  linop_handles([nX*nY, nCells+nBlobs], Af, At, 'R2R');           % Create a linear operator handle for TFOCS
    end
    
    % NOTE: In order to NOT use TFOCS, uncomment these lines and the switch
    % below (b/w TFOCS & quadprog)
%     miniBlobsMat_all{cc} = makeBlobsMatrix(nY,nX,blurRad,clipHeight);      % Create the matrix of blobs
%     H_all{cc} = 2*([rois*diag(1./prof_norms_all{cc}), miniBlobsMat_all{cc}].'*[rois*diag(1./prof_norms_all{cc}), miniBlobsMat_all{cc}]);
%     H_all{cc} = sparse(0.5*(H_all{cc}+H_all{cc}'));
%     h_all{cc} = -2*[rois*diag(1./prof_norms_all{cc}), miniBlobsMat_all{cc}].';
    
    
    switch params.dyn_mod
        case 'rwl1df'
            dp.xis    = dp.beta./lambdas_all{cc}(lambdas_all{cc}>0);
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

    dp_all{cc} = dp;
    
end



% identify blob fit to each ROI

apply_blobs = @(z) reshape(convn(reshape(z, se.movY,se.movX), one_blob, 'same'),se.movY*se.movX, 1);

% make matrix that mimics the effect of applying blobs
blobMat = zeros(size(blobsMatrix));
for bb = 1:size(blobMat,2)
    theSeed = zeros(size(blobMat,1),1);
    theSeed(bb) = 1;
    blobMat(:,bb) = apply_blobs(theSeed);
end

% get bloB fit to each cell
blobFits = nan(se.movX*se.movY,size(rois,2));
for rr = 1:size(rois,2)
    blobFits(:,rr) = lsqnonneg(blobMat,reshape(se.profiles(:,:,rr),[],1));
end



% put some parameters into sub-structs
dp = struct;
switch params.dyn_mod
    case 'rwl1df'
        dp.n_rw_iter = params.n_rw_iter;                                % Number of RWL1 iterations to run
        dp.tau       = params.tau;                                      % -----
        dp.eta       = params.eta;                                      % Parameters for RWL1: Eventually set these automatically (maybe via x-validation?)
        dp.beta      = params.beta;                                     % -----
    case 'bpdndf'
        dp.sigma2    = params.sigma2;
        dp.kappa     = params.kappa;
        dp.kappaProf = params.kappaProf;
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run per-frame analysis
    


% pick on which frames to show update
initialFrames = [1 3 10 30 100 300 1000];
initialFrames = initialFrames(initialFrames<se.movF*.2);
prcFrames     = round(linspace(0,se.movF,21));
prcFrames     = prcFrames(prcFrames>initialFrames(end));
updateFrames  = [initialFrames prcFrames];

if ~params.verbose, updateFrames = []; end


% initialize holder of frames
frameSlidingWindow = nan(se.movY*se.movX,dsT);


startTime = now;
for ff=1:nFrames
    
    % get current frame
    currFrame = reshape(se.getFrame(ff),[],1);
    
    if ff <= dsT
        % the first few times, build up a matrix of the initial frames
        frameSlidingWindow(:,ff) = currFrame;
    else
        % subequent times, delete the oldest frame and appent the current frame
        frameSlidingWindow = [frameSlidingWindow(:,2:end) currFrame];
    end
    % average the last few frames
    thisFrameFull = nanmean(frameSlidingWindow,2);
    
    
    
    % do SEUDO on each cell
    for cc = 1:length(whichCells)
        
        

        % get cell-specific info
        rois   = rois_all{cc};                                             % 
        nY     = nY_all{cc};                                               % 
        nX     = nX_all{cc};                                               % 
        nCells = sum(whichProfiles(:,cc));                                 % 
        dp     = dp_all{cc};                                               % 
        
        
        thisFrame = thisFrameFull(whichPixels(:,cc));
        
        
        % least squares solution
        tcLSQ_thisFrame = ((rois'*rois)\(rois'*thisFrame))';
        
        
        % perform SEUDO on this frame
        
        fitWeights = zeros(nCells+nY*nX,1);                                    % Initialize RWL1-DF estimate to zero
        if (ff > 1)&&(strcmp(params.dyn_mod,'bpdndf')||strcmp(params.dyn_mod,'rwl1df'))
            tc_past = [tcSEUDO(ff-1,:).';vec(fitBOB)];
            dp.k2 = update_addconst(vec(fitBOB),params,lambdas0_all{cc});
        else
            tc_past = [];
        end
        switch params.dyn_mod
            case 'rwl1df'
                fitWeights = singleRWL1DFStep(thisFrame, tc_past, ...
                    fitWeights, one_blob, rois, lambdas_all{cc}, ...
                    prof_norms_all{cc}, [nX,nY,nCells], ls_opt_all{cc}, dp); % Update weights using RWL1DF
            case 'bpdndf'
                fitWeights = singleBPDNDFStep(thisFrame, tc_past, ...
                    fitWeights, one_blob, rois, lambdas_all{cc}, dp, ...
                    prof_norms_all{cc}, [nX,nY,nCells], ls_opt_all{cc});
            otherwise
                fitWeights = solver_L1RLS(Xmat_all{cc}, thisFrame, lambdas_all{cc}, fitWeights , ls_opt_all{cc});    % get best fit weights
                % NOTE: In order to NOT use TFOCS, uncomment this line
                % and block above to make H_all etc.
                % fitWeights = quadprog(H_all{cc}, (h_all{cc}*thisFrame+lambdas_all{cc}), [], [], [], [], zeros(nX*nY+nCells,1),...
                %Inf*ones(nX*nY+nCells,1), fitWeights, optimset('Display', 'off','TolFun',1e-14,'Algorithm','trust-region-reflective')); 

        end
        fitWeights = fitWeights./normFactors_all{cc}';                                 % undo the norm factors that were needed to condition A
        
        % separate weight vector into found cell activation values (X) and
        % not-found cell activation values (Bunch Of Blobs)
        fitX   = fitWeights(1:nCells);
        fitBOB = fitWeights(nCells+1:end);
        
        % decide whether this frame is best fit by LSQ or the cells plus blobs
        
        % compute cost of each fit
        lsqCost = norm(rois*tcLSQ_thisFrame'-thisFrame,2)^2;
        if (ff > 1)&&(strcmp(params.dyn_mod,'bpdndf')||strcmp(params.dyn_mod,'rwl1df'))
            bobCost = dynamicSEUDOcost(fitWeights, thisFrame, tc_past, rois, lambdas_all{cc}, dp, params, apply_blobs_all{cc});
        else
            % Cost w/ blobs is a log-likelihood - in this case the LASSO cost [ k2 = 2*sigma2*(log(1./q-1) + sum(log(lambdas0((lambdas0>0))))) ]
            bobCost = norm(thisFrame-rois*fitX-apply_blobs_all{cc}(fitBOB),2)^2 ...
                + sum(abs(dp.k1(:).*fitWeights(:)))-dp.k2;
        end
        
        % choose the fit with lower cost
        if lsqCost<bobCost                                                 % IF LEAST-SQUARES HAS A LOWER COST
            fitFancy = tcLSQ_thisFrame;                                    % Save the least-squares solution
            fitBOB   = 0*fitBOB;                                           % Set the blob costs back to zero
            fitX     = 0*fitX;                                             % Set the blob costs back to zero
        else                                                               % IF THE SEUDO MODEL HAS A LOWER COST
            fitFancy = fitX;                                               % Save the SEUDO model solution
        end
        
        
        
        % store fits in time course variables
        tcSEUDO(ff,cc)          = fitFancy(cellIndexWithinTheseCells{cc});
        tcCellsWithBlobs(ff,cc) = fitX(cellIndexWithinTheseCells{cc});
        tcLSQ(ff,cc)            = tcLSQ_thisFrame(cellIndexWithinTheseCells{cc});
        
        if params.saveOtherCellTimeCourses
            tcOtherCellsLSQ_all{cc}(ff,:) = tcLSQ_thisFrame;
        end
        
        % save blobs
        if params.saveBlobTimeCourse
            tcBlobs{cc}(ff,:) = fitBOB;
        else
            tcBlobs(ff,1,cc)  = sum(fitBOB);
            prof              = se.profiles(:,:,cc) > 0;
            tcBlobs(ff,2,cc)  = sum(fitBOB(prof(whichPixels(:,cc))));
        end
        
        % store the costs
        lsqCosts(ff,cc) = lsqCost;
        bobCosts(ff,cc) = bobCost;
        
        
        
        
        % store updated dp
        dp_all{cc} = dp;
        
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
extras.lsqCosts         = lsqCosts;
extras.blobCosts        = bobCosts;
extras.tcCellsWithBlobs = tcCellsWithBlobs;
extras.tcBlobs          = tcBlobs;
extras.one_blob         = one_blob;
extras.blobFits         = blobFits;
extras.whichPixels      = whichPixels;
extras.whichProfiles    = whichProfiles;
extras.tcOtherCellsLSQ  = tcOtherCellsLSQ_all;












% STORE RESULTS IN SEUDO OBJECT



idx = length(se.tcSeudo) + 1;                           % identify index to store in


se.tcSeudo(idx).tc = tcSEUDO;                           % SEUDO time courses 
se.tcSeudo(idx).tcLSQ = tcLSQ;                          % LSQ time courses 
se.tcSeudo(idx).params = params;                        % parameters
se.tcSeudo(idx).extras = extras;                        % intermediate computations for debugging



if params.verbose
    fprintf('   completed in %0.1f min (%s), results stored in field tcSeudo(%d)\n',...
        (now-startTime)*24*60,datestr(now),idx)
end



end












% SUBFUNCTIONS



function A = makeBlobsMatrix(nY,nX,blurRadius,clipHeight)


% make a single blob
fblur = fspecial('gauss',[2*nY+1 2*nX+1],blurRadius);

% zero any pixels below a certain height
fblur = fblur .* double(fblur > clipHeight * max(fblur(:)));

% make copy of blob centered at each pixel
A = nan(nX*nY);
for xx=1:nX
    for yy=1:nY
        A(sub2ind([nY nX],yy,xx),:) = reshape(fblur(nY-yy+1+(1:nY),nX-xx+1+(1:nX)),1,[]);
    end
end

% give each blob the same sum
A = bsxfun(@rdivide,A,sum(A,2));


end





function bobCost = dynamicSEUDOcost(fitWeights, dF_now, tc_old, rois, lambdas, dp, params, apply_blobs)

% bobCost = dynamicSEUDOcost(fitWeights, dF_now, tc_old, rois, lambdas, dp, p, apply_blobs)
%
%
%
% 2018 - Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some initializations

nRois = numel(fitWeights)-numel(dF_now);

fitX   = fitWeights(1:nRois);
fitBOB = fitWeights(nRois+1:end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate cost

switch params.dyn_mod
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
%%%%%%%%%%%%%%%%%%





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




