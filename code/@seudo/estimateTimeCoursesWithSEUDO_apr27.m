function se = estimateTimeCoursesWithSEUDO(se,varargin)
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
%   dsTime - # frames to average. SEUDO analysis is applied to the average of the most recent <dsTime> frames 
%   lambda - sparsity multiplicative value
%   lambdaProf - profile sparsity multiplicative value
%   lambdaBlob - blob sparsity multiplicative value
%   blobRadius - radius (in pixels) of a single Gaussian blob
%
%   whichCells - list of which cells to perform analysis on. if empty (default), use all.
%                   unselected cells have a time course of NaNs.
%
%   verbose - logical, whether to print progress to the command window
%   showDefaults - if true, display default values and return without performing SEUDO
%   validateParams - if true, validate parameters without performing SEUDO
%                       if all parameters are valid, returns empty
%                       if any parameters are invalid, returns the first error
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
%               tcBlobs - (# frames) x (# pixels) matrix, each column is the time course of one blob
%               ...
%   
%


% Get Input Parameters

% parase parameters
p = inputParser;
p.addParameter('p'          , 0.00001, @(x)x>=0&&x<=1);   
p.addParameter('sigma2'     , .01);   
p.addParameter('dsTime'     , 1);
% p.addParameter('lambda'     , 1/50);
p.addParameter('lambdaProf' , 0);
p.addParameter('lambdaBlob' , 20); 
p.addParameter('blobRadius' , 1.2);

p.addParameter('dyn_mod'    , 'static');

p.addParameter('kappa'      , 0.05);
p.addParameter('kappaProf'  , 0.001);
p.addParameter('tau'        , 0.1);
p.addParameter('eta'        , 0.01);
p.addParameter('beta'       , 0.01);
p.addParameter('n_rw_iter'  , 3);

p.addParameter('verbose'        , true);
p.addParameter('validateParams' , false);
p.addParameter('showDefaults'   , false);
p.addParameter('whichCells' , []);



% just validate inputs, if desired
q               = inputParser;
q.KeepUnmatched = true;
q.addParameter('validateParams' , false, @(x)isscalar(x)&islogical(x));
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
q.addParameter('showDefaults' , false, @(x)isscalar(x)&islogical(x));
parse(q,varargin{:});
if q.Results.showDefaults
    % display default values
    fprintf('\ndefault values for SEUDO:\n\n')
    parse(p);
    se = p.Results;
    disp(se)
    %fprintf('\n')
    return
end


% parse inputs
parse(p,varargin{:});


% put some parameters into sub-structs
switch p.Results.dyn_mod
    case 'rwl1df'
        dp.n_rw_iter = p.Results.n_rw_iter;                                % Number of RWL1 iterations to run
        dp.tau       = p.Results.tau;                                      % -----
        dp.eta       = p.Results.eta;                                      % Parameters for RWL1: Eventually set these automatically (maybe via x-validation?)
        dp.beta      = p.Results.beta;                                     % -----
    case 'bpdndf'
        dp.sigma2    = p.Results.sigma2;
        dp.kappa     = p.Results.kappa;
        dp.kappaProf = p.Results.kappaProf;
end



% copy to struct that will be output
params = p.Results;

% set code version
params.version = 1;

% ESTIMATE TIME COURSES 

if isempty(params.whichCells)
    whichCells = 1:size(se.profiles,3);
else
    whichCells = params.whichCells;
end


% Get Useful Sizes

% note original movie size
nY = se.movY;                                                              % Get size of movie
nX = se.movX;                                                              % -----------------

% note sizes
nFrames = se.movF;                                                         % Number of time-frames
nRois   = length(whichCells);                                              % Number of ROIs
nBlobs  = nX*nY;                                                           % Number of blobs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set some parameters


dsT    = params.dsTime;                                                    % temporal downsample
p      = params.p;                                                         % p = prob(not found cell is active)
sigma2 = params.sigma2/dsT;                                                % noise variance

% reshape input data into matrices
rois = reshape(se.profiles(:,:,whichCells),[],nRois);                      % nPix x nRois
%dF = reshape(dF,[],size(dF,3));                                           % nPix x nFrames


lambdas0 = [params.lambdaProf*ones(1,nRois), ...
                                      params.lambdaBlob*ones(1,nBlobs)]';  % set lambda values

% set constants automatically
k1 = 2*sigma2*lambdas0;
k2 = 2*sigma2*(log(1./p-1) + sum(log(lambdas0(lambdas0>0))));
k3 = 2*sigma2*(log(1./p-1));

% initialize storage of results
tcLSQ            = nan(nFrames,nRois);
tcSEUDO          = nan(nFrames,nRois);                                     % final output of the SEUDO-LS comparison
lsqCosts         = nan(nFrames,1);                                         % Least-squares cost
bobCosts         = nan(nFrames,1);                                         % Total cost of SEUDO cell construction
tcCellsWithBlobs = nan(nFrames,nRois);                                     % Total cost for the known cells and the SEUDO cell construction 
tcBlobs          = nan(nBlobs, nFrames);                                   % Total cost of SEUDO cell construction
fitGraphic       = [];                                                     % 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set solver options
ls_opt            = struct;                                                % Create TFOCS options struct
ls_opt.tol        = 0.01;                                                  % Set algorithmic tolerance (for TFOCS)
ls_opt.nonneg     = true;                                                  % Set non-negative flag (for TFOCS)
ls_opt.printEvery = 0;                                                     % Suppress TFOCS output


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make one blob

% parameters
blurRad    = params.blobRadius;
clipHeight = 0.0;

blobsMatrix = makeBlobsMatrix(nY,nX,blurRad,clipHeight);                   % Create the matrix of blobs
one_blob    = reshape(blobsMatrix(:,...                                    
                   sub2ind([nY nX],round(0.5*nY),round(0.5*nX))),nY,nX).'; % Isolate a single blob - start with a single frame with the blob in the middle ...
[imax,jmax] = find(one_blob==max(one_blob(:)),1);                          % ... Find the blob center ...
one_blob    = one_blob(imax+(-7:7),jmax+(-7:7));                           % ... And crop the blob
one_blob    = one_blob/sqrt(sum(one_blob(:).^2));
  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create TFOCS Linear Operator

prof_norms  = sqrt(sum(rois.^2,1));                                        % Get only the profile norms
normFactors = [prof_norms,ones(1,nBlobs)];
lambdas     = k1./normFactors';                                            % scale lambdas to match
apply_blobs = @(z) reshape(convn(reshape(z, nY,nX), one_blob, 'same'),...
                                                                nX*nY, 1); % Make an anonymous function that applies the blobs accross the image via convolution

if strcmp(params.dyn_mod,'static')
    Af   = @(z) (rois*diag(1./prof_norms))*vec(z(1:nRois)) + reshape(convn(...
             reshape(z(nRois+(1:nX*nY)), nY,nX), one_blob, 'same'), nX*nY, 1); % Make forward operator
    At   = @(z) [(rois*diag(1./prof_norms)).'*vec(z); reshape(convn(...
                             reshape(z, nY,nX), one_blob, 'same'), nBlobs,1)]; % Make transpose operator
    Xmat =  linop_handles([nX*nY, nRois+nBlobs], Af, At, 'R2R');               % Create a linear operator handle for TFOCS
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run per-frame analysis

switch params.dyn_mod
    case 'rwl1df'
        dp.xis    = dp.beta./lambdas(lambdas>0);
        dp.alphas = 0.5*dp.tau./dp.xis - 1;
        dp.k1     = k1;
        dp.k2     = k2;
        dp.k3     = k3;
    case 'bpdndf'
        dp.k1     = k1;
        dp.k2     = k2;
        dp.k3     = k3;
end


startTime = now;

% pick on which frames to show update
initialFrames = [1 3 10 30 100 300 1000];
initialFrames = initialFrames(initialFrames<se.movF*.2);
prcFrames     = round(linspace(0,se.movF,21));
prcFrames     = prcFrames(prcFrames>initialFrames(end));
updateFrames  = [initialFrames prcFrames];

if ~params.verbose, updateFrames = []; end


% initialize holder of frames
frameSlidingWindow = nan(se.movY*se.movX,dsT);

    
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
    thisFrame = nanmean(frameSlidingWindow,2);
    
    
    
    % least squares solution
    tcLSQ(ff,:) = ((rois'*rois)\(rois'*thisFrame))';  
    
    % get previous time course if needed
    if (ff > 1)&&(strcmp(params.dyn_mod,'bpdndf')||strcmp(params.dyn_mod,'rwl1df'))
        tc_past = [tcSEUDO(ff-1,:).';vec(fitBOB)];
    else
        tc_past = [];
    end
    
    
    % perform SEUDO on this frame
    
    fitWeights = zeros(nRois+nBlobs,1);                                    % Initialize RWL1-DF estimate to zero
    if (ff > 1)&&(strcmp(params.dyn_mod,'bpdndf')||strcmp(params.dyn_mod,'rwl1df'))
        tc_past = [tcSEUDO(ff-1,:).';vec(fitBOB)];
        k2 = update_addconst(vec(fitBOB),params,lambdas0);
    else
        tc_past = [];
    end
    switch params.dyn_mod
        case 'rwl1df'
            fitWeights = singleRWL1DFStep(thisFrame, tc_past, ...
                                 fitWeights, one_blob, rois, lambdas, ...
                                   prof_norms, [nX,nY,nRois], ls_opt, dp); % Update weights using RWL1DF
        case 'bpdndf'
            fitWeights = singleBPDNDFStep(thisFrame, tc_past, ...
                             fitWeights, one_blob, rois, lambdas, dp, ...
                                       prof_norms, [nX,nY,nRois], ls_opt);
        otherwise
            fitWeights = solver_L1RLS(Xmat, thisFrame, lambdas, fitWeights , ls_opt);    % get best fit weights
    end
    fitWeights = fitWeights./normFactors';                                 % undo the norm factors that were needed to condition A
    
    % separate weight vector into found cell activation values (X) and
    % not-found cell activation values (Bunch Of Blobs)
    fitX   = fitWeights(1:nRois);
    fitBOB = fitWeights(nRois+1:end);
    
    % decide whether this frame is best fit by LSQ or the cells plus blobs
    
    % compute cost of each fit
    lsqCost = norm(rois*tcLSQ(ff,:)'-thisFrame,2)^2;
    if (ff > 1)&&(strcmp(params.dyn_mod,'bpdndf')||strcmp(params.dyn_mod,'rwl1df'))
        bobCost = dynamicSEUDOcost(fitWeights, thisFrame, tc_past, rois, lambdas, dp, params, apply_blobs);
    else
        bobCost = norm(thisFrame-rois*fitX-apply_blobs(fitBOB),2)^2 ...
                                      + sum(abs(k1(:).*fitWeights(:)))-k2; % Cost w/ blobs is a log-likelihood - in this case the LASSO cost [ k2 = 2*sigma2*(log(1./q-1) + sum(log(lambdas0((lambdas0>0))))) ]
    end
    
    % choose the fit with lower cost
    if lsqCost<bobCost                                                     % IF LEAST-SQUARES HAS A LOWER COST
        fitFancy = tcLSQ(ff,:);                                            % Save the least-squares solution
        fitBOB   = 0*fitBOB;                                               % Set the blob costs back to zero
        fitX     = 0*fitX;                                                 % Set the blob costs back to zero
    else                                                                   % IF THE SEUDO MODEL HAS A LOWER COST
        fitFancy = fitX;                                                   % Save the SEUDO model solution
    end
    
    
    % store fits in time course variables
    tcSEUDO(ff,:)          = fitFancy;
    tcCellsWithBlobs(ff,:) = fitX;
    tcBlobs(:,ff)          = fitBOB;
    
    % store the costs
    lsqCosts(ff) = lsqCost;
    bobCosts(ff) = bobCost;
    
    
    
    %if mod(ff,100)==0
    if ismember(ff,updateFrames)
        finishDur = (nFrames-ff)/ff*(now-startTime);
        %fprintf('Time-step %d of %d finished in %f s.\n', ff, nFrames, tframe)
        fprintf('      completed %d of %d (%0.1f%%) in %0.1f min, estimated finish in %0.1f min (%s)\n',...
            ff, nFrames, 100*ff/nFrames,(now-startTime)*24*60,finishDur*24*60,datestr(now+finishDur))
    end
    
    
    
end

extras = struct;
extras.lsqCosts = lsqCosts;
extras.blobCosts = bobCosts;
extras.tcCellsWithBlobs = tcCellsWithBlobs;
extras.fancyFits = fitGraphic;
extras.tcBlobs = tcBlobs;
extras.one_blob = one_blob;













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
                                      + sum(abs(k1(:).*fitWeights(:)))-k2; % Cost w/ blobs and dynamics is ...
    otherwise
        error('dynamics type %s not recoanized',dyn_mod)
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%





function k2 = update_addconst(fitBoB, p, lambdas0)

% k2 = update_addconst(fitBoB, p, lambdas0)
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
        k2 = 2*params.sigma2*(log(1./params.p-1) + sum(log(lambdas0((lambdas0>0)))));
    case 'bpdndf'  
        logZ = params.lambdaBlob*sum(abs(fitBoB)) ...
            - params.lambdaBlob*numel(fitBoB)./(4*params.kappa) ...
            - sum(log(normcdf(sqrt(2*params.kappa)*fitBoB - ...
                    params.lambdaBlob/sqrt(2*params.kappa))));
        k2 = 2 *params.sigma2*(log(1./params.p-1) - logZ );
    otherwise
        k2 = 2*params.sigma2*(log(1./params.p-1) + sum(log(lambdas0((lambdas0>0)))));
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
        Xmat =  linop_handles([nX*nY, nRois+nBlobs], Af, At, 'R2R');       % Create a linear operator handle for TFOCS

        fitWeights = solver_L1RLS(Xmat,dF_now,lambdas,...
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
    Xmat  = linop_handles([nX*nY, nRois+nBlobs], Af, At, 'R2R');           % Create a linear operator handle for TFOCS
    fitWeights = solver_L1RLS(Xmat,dF_now,lambdas,...
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




