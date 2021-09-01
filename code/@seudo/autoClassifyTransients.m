function varargout = autoClassifyTransients(se,structSpec,varargin)
% seudoObject.autoClassifyTransients(structSpec,...)
%
% automatically classify transients
%
%
%
% EXAMPLES:
%
%   se.autoClassifyTransients('default')
%   se.autoClassifyTransients('default','overwrite',true,'method','lbq')
%   se.autoClassifyTransients({'seudo',3},'radius',[2 6])
%
%   [transientMetricsByCell,params] = seudoObject.autoClassifyTransients(...)
%
%
% INPUTS:
%
%   structSpec - which time course struct to classify, e.g. 'default', {'seudo',3}
%           for options, see help for pickTcStruct
%
%
% PARAMETERS:
%
%   whichCells - which cells to classify, if empty use all cells
%   saveResults - if true (default), store results in the seudo object
%                 if false, return results in a variable
%   overwrite - overwrite existing classification, only applies when saveResults is true
%                   if false (default), will return error if any cells are already classified
%   plotComparison - if true (default = false), compare existing classification to new automated classification
%                       only applies if saveResults is false
%   profileFieldName - name of field in seudoObject.tcDefault.transientInfo to get transient profiles, default = 'shapes'
%   method - which classification method to use
%               'lbq' - see paper
%               'corr' - classify as true if Pearson correlation between transient profile and source profile exceeds threshold
%   
% 
%   PARAMETERS FOR METHOD 'lbq':
%
%           radius - 2-length vector, inner and outer radius of autocorrelation lags to include, units = pixels
%           tauBounds - 2-length vector, bounds for marking a transient as mixed, see text (tau_1 and tau_2)
%           ignoreZeros - ignore any pixels that have zero profile weight (default = false)
%
%
%   PARAMETERS FOR METHOD 'corr':
%
%           corrThresh - threshold value for classifying transient true
%
%
%
%
% OUTPUTS:
%
%   if saveResults is true,  results are stored in these fields of the time course struct:
%
%       transientInfo.classification - standard classification structure (see help for computeTransientInfo)
%       transientInfo.autoClass - struct of fields with intermediate results of the classification (see below)
%
%
%   if save results is false, output is like this:
%       [transientMetricsByCell,params] = seudoObject.autoClassifyTransients(...)
%
%   transientMetricsByCell - struct array, length equal to the number of analyzed cells
%
%       for method 'lbq', fields are:
%
%           acm - T-length vector, ????
%           lbq - T-length binary vector, how each transient was classified (0 = not contaminated, 1 = contaminated)
%           res - Y x X x T matrix, residual after subtracting the best fit of the cell profile from each transient shape
%           corrs - Y x X x T matrix, autocorrelation of each residual
%           lbq_ext - ?????
%
%           for T = # transients, X = # pixels in X, Y = # pixels in Y
%
%
%       for method 'corr', fields are:
%
%           acm - T-length vector, ????
%
%
%   params - struct of parameters
%
%
% 2018 Jeff Gauthier & Adam Charles





% parse inputs
p = inputParser;
p.addParameter('whichCells',[]);
p.addParameter('overwrite',false);
p.addParameter('saveResults',true);
p.addParameter('plotComparison',false);

p.addParameter('profileFieldName','shapes');
p.addParameter('ignoreZeros',false);
p.addParameter('minimumWindow',true);

p.addParameter('blurRadius',1);
p.addParameter('blurProfilesForFitting',false);

p.addParameter('method','corr');
p.addParameter('corrThresh',0.4);

p.addParameter('radius',[4 6]);
p.addParameter('alpha',0.01);
p.addParameter('dof',0);
p.addParameter('tauBounds',[0.2,0.25]);

p.addParameter('weightedTrans',false);
p.addParameter('resRatioMean',true);
parse(p,varargin{:});


params = p.Results;
clear p

% parse user input
structSpec = pickTcStruct(se,structSpec);


% ensure transient info is computed
if ~isfield(se.(structSpec{1})(structSpec{2}),'transientInfo') || isempty(se.(structSpec{1})(structSpec{2}).transientInfo)
    startTime = now;
    fprintf('computing transient info...')
    se.computeTransientInfo(structSpec);
    fprintf('done in %0.1f min\n',(now-startTime)*24*60)
end

% parameters

r1 = params.radius(1);
r2 = params.radius(1+mod(length(params.radius)-1,2));

if isempty(params.whichCells)
    whichCells = 1:se.nCells;
else
    whichCells = params.whichCells;
end


% confirm before overwriting

if params.saveResults
    % identify which cells already have classification
    hasClass = cellfun(@(x)~all(isnan(x)),{se.(structSpec{1})(structSpec{2}).transientInfo.classification});
    
    if any(hasClass(whichCells)) && ~params.overwrite
        error('classification exists for %d of %d cells, to overwrite call with arguments: autoClassifyTransients(...,''overwrite'',true)',...
            sum(hasClass(whichCells)),length(whichCells))
    end
end


for cc = 1:length(whichCells)
    
    % get info for this cell
    
    cellID     = whichCells(cc);                                           % Get the current cell ID
    
    metric_out = struct;                                                   % initialize output struct
    
    T = se.(structSpec{1})(structSpec{2}).transientInfo(cellID);           % Pull out the current transient information for this cell
    
    if size(T.times,1) == 0                                                % If there are no transients in this time range, return appropriate struct
        
        
        metric_out.cfrac  = nan;                                           % primary metric that counts the fraction of transients that are false
        metric_out.rfrac  = nan;                                           % secondary metric that counts the relative weight of false transient residuals  
        
        switch params.method
            case 'lbq'
                metric_out.acm    = [];                                    % Initialize output array of accumulated correlation values
                metric_out.lbq    = [];                                    % Initialize output array of Ljung-Box test results
                metric_out.relres = [];                                    % Initialize output array of residual norms
                metric_out.res    = [];                                    % Initialize the array of residuals
                metric_out.corrs  = [];                                    % Initialize the array of correlation coefficients
                metric_out.classification = [];
            case 'corr'
                metric_out.corrs = [];
            otherwise
                error('classification method ''%s'' not recognized',params.type)
        end
            
        % store in seudo object
        if params.saveResults
            
            % store correlation
            se.(structSpec{1})(structSpec{2}).transientInfo(cellID).classification = [];

            % store metric (with parameters)
            metric_out.params = params;
            se.(structSpec{1})(structSpec{2}).transientInfo(cellID).autoClass = metric_out;

        end
        
        % store in variable to be returned
        if nargout > 0
            
            if cc == 1
                transientMetricsByCell = metric_out;
            else

                if isstruct(transientMetricsByCell)
                    [metric_out,transientMetricsByCell] = matchOtherFields(metric_out,transientMetricsByCell);
                end
                transientMetricsByCell(cc) = metric_out;

            end

        end
        
        continue
    end
    
    tc      = se.(structSpec{1})(structSpec{2}).tc(:,cellID);              % get full time course
    evts    = T.times;                                                     % transient time spans for significant transients
    tc_amps = zeros(size(evts,1),1);                                       % Initialize array of amplitudes
    for ll = 1:size(evts,1)
        tc_amps(ll) = max(tc(evts(ll,1):evts(ll,2)));
    end
    x_prof  = se.profiles(T.window(1):T.window(2),...                      % get profile for this cell
                                         T.window(3):T.window(4), cellID);
    x_other = se.profiles(T.window(1):T.window(2),...                      % also get the profiles for other cells (but keep them separate)
                       T.window(3):T.window(4),[1:cellID-1,cellID+1:end]);
    
    
    if params.minimumWindow
        % cut to the smallest bounding rectangle
        winSubY = any(x_prof>0,2);
        winSubX = any(x_prof>0,1);
    else
        % use all pixels in the shape calculation
        winSubY = true(size(x_prof,1),1);
        winSubX = true(size(x_prof,2),1);
    end
    
    % keep only the pixels in the window
    x_prof  = x_prof(winSubY,winSubX);                                     % Restrict the cell profile to the given window...
    Ffits   = T.(params.profileFieldName)(winSubY,winSubX,:);              %    ... as well as the data...
    x_other = x_other(winSubY,winSubX,:);                                  %    ... as well as the other cell profiles
    
    
    % note size of final window
    wY = sum(winSubY);                                                     % 
    wX = sum(winSubX);                                                     % 
    
    
    if params.blurProfilesForFitting
        % blur profiles before fitting
        x_profForFit    = imgaussfilt(x_prof,params.blurRadius);
        x_otherForFit   = imgaussfilt(x_other,params.blurRadius);
    else
        % don't blur profiles before fitting
        x_profForFit = x_prof;
        x_otherForFit = x_other;
    end
    
    x_profForFit    = x_profForFit(:);
    x_otherForFit   = reshape(x_otherForFit,wX*wY,[]);
    x_otherForFit = x_otherForFit(:,sum(x_otherForFit)~=0);
    
    
    % reshape
    x_prof  = x_prof(:);
    Ffits_nofilt = reshape(Ffits,wX*wY,[]);
    x_other = reshape(x_other,wX*wY,[]);
    x_other = x_other(:,sum(x_other)~=0);
    
    
    % compute sum of all transients
    tot_tcs = 0;
    for kk = 1:size(evts,1)
        tot_tcs = tot_tcs + sum(tc(evts(kk,1):evts(kk,2)));
    end
    
    
    % Low-pass filtering
    Ffits        = imgaussfilt(Ffits,params.blurRadius);                   % Low-pass filter Fluorescence data with a Gaussian filter
    Ffits        = reshape(Ffits,wX*wY,[]);                                % Reshape Fluorescence array to be a matrix with each frame as a column
    
    
    % Calculate fit
    fit_vals         = [x_profForFit,x_otherForFit]\Ffits;                 % Calculate the fits to all ROIs in the image
    
    %%%% Block of code that can implement non-negative fits -- ASC %%%%%%%%
%     XA = [x_prof,x_other];
%     fit_vals = zeros(size(XA,2), size(Ffits,2));
%     for ll = 1:size(Ffits,2)
%         fit_vals = quadprog(2*(XA'*XA), -2*XA'*Ffits(:,ll), [], [], [], [], ...
%                           zeros(size(XA,2),1),[],[],optimset('Display','off'));
%     end
    %%%% END :: Block of code that can implement non-negative fits %%%%%%%%
    
    res_evt          = Ffits - x_otherForFit*fit_vals(2:end,:);            % Calculate event-triggering profile: residual images removing all but ROI in question
    
    res              = Ffits - [x_profForFit,x_otherForFit]*fit_vals;      % Calculate the residual images - based on all ROIs
    res              = reshape(res,[wY,wX,size(evts,1)]);                  % Reshape residuals
    
    metric_out.resNorms  = reshape(sum(sum(res.^2,1),2),[],1);
    metric_out.resNoCell = reshape(sum(res_evt.^2,1),[],1);
    metric_out.resRatios = reshape(sum(sum(res.^2,1),2),[],1)./reshape(sum(res_evt.^2,1),[],1);
    metric_out.residuals = res;
    metric_out.residualsWithoutSource = reshape(res_evt,[wY,wX,size(evts,1)]);
    
    Ffits = Ffits_nofilt;
    
    metric_out.cfrac  = 0;                                                 % primary metric that counts the fraction of transients that are false
    metric_out.rfrac  = 0;                                                 % secondary metric that counts the relative weight of false transient residuals  
    metric_out.rfrac2 = 0;
    metric_out.rfrac3 = 0;
    
    
    res_off  = 0;
    
    switch params.method
        
        case 'lbq'
            
            if params.ignoreZeros                                          % set to nan any pixels that have zero profile weight
                zeroPix = x_prof==0;
            else
                zeroPix = [];
            end
            
            [gX, gY]          = meshgrid(-(wX-1):(wX-1), -(wY-1):(wY-1));  % Create a meshgrid of points to isolate the required correlations
            metric_out.acm    = zeros(size(evts,1),1);                     % Initialize output array of accumulated correlation values
            metric_out.lbq    = zeros(size(evts,1),1);                     % Initialize output array of Ljung-Box test results
            metric_out.lbq_ext = zeros(size(evts,1),3);
            metric_out.relres = zeros(size(evts,1),1);                     % Initialize output array of residual norms
            metric_out.res    = zeros(wY,wX, size(evts,1));                % Initialize the array of residuals
            metric_out.corrs  = zeros(size(gX,1),size(gX,2),size(evts,1)); % Initialize the array of correlation coefficients
            
%             res_den           = sum(nans2zeros(res_evt(:)).^2);
            res_den           = sum(nans2zeros(res_evt).^2,1);
            for kk = 1:size(evts,1)                                        % Get auto-correlation for each variable
                TMPRES          = res(:,:,kk);                             % Temporarily store the residual (in case of requiring modification)
                if params.ignoreZeros
                    TMPRES(zeroPix) = nan;                                 % set to nan any pixels that have zero profile weight
                end
                c               = xcorr2(nans2zeros(TMPRES));              % Calculate the autocorrelation values
                c               = c/abs(c((gX==0)&(gY==0)));               % Normalize the autocorrelation to get coefficients
                metric_out.corrs(:,:,kk) = c;                              % Store the autocorrelation values
                metric_out.res(:,:,kk)   = TMPRES;                         % Store the residual
                
                c(sqrt(gX.^2 + gY.^2) < r1) = 0;                           % Remove inner correlation coefficients
                c(sqrt(gX.^2 + gY.^2) > r2) = 0;                           % Remove outer correlation coefficients
                metric_out.acm(kk)          = norm(c(:));                  % Save the norm of the annulus of correlation coefficients
                
                if isempty(params.dof)
                    dof_val = ceil(numel(x_prof>0));
                else
                    dof_val = params.dof;
                    if dof_val < 1
                        dof_val = ceil(dof_val*numel(x_prof>0));
                    end
                end
                
                [metric_out.lbq(kk), metric_out.lbq_ext(kk,1), ...         % Calculate the 2D Ljung-Box test
                metric_out.lbq_ext(kk,2), metric_out.lbq_ext(kk,3)] = ...
                    lbqtest2d(double(TMPRES),...                           
                    'Alpha',params.alpha,'Lags',[r1,r2], 'DoF', dof_val);  % Ljung-Box test parameters
                
                metric_out.relres(kk) = (norm(nans2zeros(TMPRES(:)))./norm(res_evt(:,kk))).^2;
                if metric_out.lbq(kk) == 1
                    if params.weightedTrans
%                         avgWeight = fit_vals(1,kk)/sum(fit_vals(1,:));     % average weight is based on transient size
                        avgWeight = tc_amps(kk)/sum(tc_amps);              % average weight is based on transient size
                    else 
                        avgWeight = 1./size(evts,1);                       % average weight is just based on the number of events total (all transients have the same weight)
                    end
                    metric_out.cfrac = metric_out.cfrac + avgWeight;       % Save the primary metric that counts the fraction of transients that fail the 2D LB Q-test
                    if params.resRatioMean
                        metric_out.rfrac = metric_out.rfrac + ...               
                                    sum(nans2zeros(TMPRES(:)).^2); % secondary metric that counts the relative weight of false transient residuals   
                    else
                        metric_out.rfrac = metric_out.rfrac + ...             
                           sum(nans2zeros(TMPRES(:)).^2)/(res_off + sum(res_evt(:,kk).^2)); % secondary metric that counts the relative weight of false transient residuals       
                    end
                end
            end
            if ~params.resRatioMean
                if sum(metric_out.lbq==1) > 0
                    metric_out.rfrac = metric_out.rfrac/(res_off + sum(metric_out.lbq==1)); 
                end
            else
                if sum(metric_out.lbq==1) > 0
                    metric_out.rfrac = metric_out.rfrac/sum(res_den(metric_out.lbq==1)); 
                end
            end
            % use the metric to create a classification
            newClassification = seudo.valTrue*ones(size(T.times,1),1);
            newClassification((metric_out.lbq(:)==1)&(metric_out.relres(:)>params.tauBounds(2))) = seudo.valFalse;
            newClassification((metric_out.lbq(:)==1)&(metric_out.relres(:)>params.tauBounds(1))&(metric_out.relres(:)<params.tauBounds(2))) = seudo.valMix;
            
            metric_out.fullVal = sum(tc_amps.*reshape(metric_out.resRatios,[],1))/sum(tc_amps);
            
        case 'corr'
            
            % compute corr for each cell
            
            if params.ignoreZeros                                          % set to nan any pixels that have zero profile weight
                zeroPix = x_prof==0;
                x_prof(zeroPix) = nan;
                Ffits(zeroPix,:) = nan;
            end
           
            corrs = nan(size(T.times,1),1);                                % initialize output
            for tt = 1:size(T.times,1)          
                corrs(tt) = corr(x_prof,Ffits(:,tt),'rows','pairwise');    % compute correlation, ignoring nans
            end 
            metric_out.corrs = corrs;                                      % store correlation values
            
            newClassification = seudo.valFalse*ones(size(T.times,1),1);    % set all transients to false
            newClassification(corrs >= params.corrThresh) = seudo.valTrue; % mark true any with sufficiently high correlation
            
            
            res_den_all  = sum(nans2zeros(res_evt).^2,1);
            res_den2 = sum(res_den_all(newClassification==seudo.valTrue));
            res_den  = sum(res_den_all(newClassification==seudo.valFalse));
            for kk = 1:size(evts,1)                                        % loop over transients
                if ismember(newClassification(kk),[seudo.valFalse seudo.valTrue])
                    switch newClassification(kk)
                        case seudo.valFalse
                            rfracFieldName = 'rfrac';
                            rd=res_den;
                        case seudo.valTrue
                            rfracFieldName = 'rfrac2';
                            rd=res_den2;
                        otherwise
                            error('code bug')
                    end
                    TMPRES          = res(:,:,kk);                         % Temporarily store the residual
                    if params.ignoreZeros
                        TMPRES(zeroPix) = nan;                             % set to nan any pixels that have zero profile weight
                    end
                    if params.weightedTrans
%                         avgWeight = fit_vals(1,kk)/sum(fit_vals(1,:));     % average weight is based on transient size
                        avgWeight = tc_amps(kk)/sum(tc_amps);              % average weight is based on transient size
                    else 
                        avgWeight = 1./size(evts,1);                       % average weight is just based on the number of events total (all transients have the same weight)
                    end
                    if ismember(newClassification(kk),[seudo.valFalse])
                        metric_out.cfrac = metric_out.cfrac + avgWeight;       % primary metric that counts the fraction of transients that are false
                    end
                    if params.resRatioMean
                        metric_out.(rfracFieldName) = metric_out.(rfracFieldName) + ...
                            sum(nans2zeros(TMPRES(:)).^2)/rd; % secondary metric that counts the relative weight of false transient residuals
                        metric_out.rfrac3 = metric_out.rfrac3  + ...
                            sum(nans2zeros(TMPRES(:)).^2)/sum(res_den_all);
                    else
                        metric_out.(rfracFieldName) = metric_out.(rfracFieldName) + ...
                            sum(nans2zeros(TMPRES(:)).^2)/sum(res_evt(:,kk).^2); % secondary metric that counts the relative weight of false transient residuals
                        metric_out.rfrac3 = metric_out.rfrac3  + ...
                            sum(nans2zeros(TMPRES(:)).^2)/sum(res_evt(:,kk).^2);
                    end
                    
                end
                
            end
            if ~params.resRatioMean
                if sum(newClassification==seudo.valFalse) > 0
                    metric_out.rfrac = metric_out.rfrac/sum(newClassification==seudo.valFalse);
                end
                if sum(newClassification==seudo.valTrue) > 0
                    metric_out.rfrac2 = metric_out.rfrac2/sum(newClassification==seudo.valTrue);
                end
                metric_out.rfrac3 = metric_out.rfrac3/numel(newClassification); 
            end
            
            metric_out.fullVal = sum(tc_amps.*reshape(metric_out.resRatios,[],1))/sum(tc_amps);
            % -------------------
            
            
        otherwise
            error('classification method ''%s'' not recognized',params.type)
            
    end
    
    metric_out.classification = newClassification;
    
    
    % save params in output
    metric_out.params = params;
        
    % store in seudo object
    if params.saveResults
        
        % store correlation
        se.(structSpec{1})(structSpec{2}).transientInfo(cellID).classification = newClassification;
        
        % store metric in seudo object
        se.(structSpec{1})(structSpec{2}).transientInfo(cellID).autoClass = metric_out;
        
    end
    
    % store in variable to be returned
    if nargout > 0
        
        if cc == 1
            transientMetricsByCell = metric_out;
        else
            if isstruct(transientMetricsByCell)
                    [metric_out,transientMetricsByCell] = matchOtherFields(metric_out,transientMetricsByCell);
            end
            % for some inexplicable reason, this is necessary
            fn = fieldnames(metric_out);
            for ff = 1:length(fn)
                transientMetricsByCell(cc).(fn{ff}) = metric_out.(fn{ff});
            end
            
            %transientMetricsByCell(cc) = metric_out; % would be nice if this worked!
        end
        
    end
    
end




% return results as variables
if nargout > 0
    varargout{1} = transientMetricsByCell;
end


% plot comparison, if desired
if ~params.saveResults && params.plotComparison
    % note: this assumes the classification in "se.tcDefault" was done by a human
    
    
    % gather classifications from all transients, all cells
    classHuman = [];
    classAuto = [];
    for cc = 1:se.nCells
        classHuman = cat(1,classHuman,se.tcDefault.transientInfo(cc).classification);
        classAuto  = cat(1,classAuto,transientMetricsByCell(cc).classification);
    end
    
    % count up number of transients in each category
    Ntrue = hist(classHuman(classAuto == seudo.valTrue),[seudo.valFalse seudo.valMix seudo.valTrue]);
    Nfalse = hist(classHuman(classAuto == seudo.valFalse),[seudo.valFalse seudo.valMix seudo.valTrue]);
    Ntrue(end+1) = sum(isnan(classHuman(classAuto == seudo.valTrue)));
    Nfalse(end+1) = sum(isnan(classHuman(classAuto == seudo.valFalse)));
    
    
    % plot comparison
    figure
    bar([Ntrue; Nfalse]',1)
    set(gca,'xtick',1:4,'xticklabel',{'false','mixed','true','unclassified'},'tickdir','out')
    ylabel('# transients')
    xlabel('existing classification')
    colormap([.5 .5 1; 1 .5 .5])
    legend('auto: true','auto: false')
    box off
    
    
    % alternative plotting
    if 0 || 0
        figure
        subplot(2,1,1)
        bar(Ntrue)
        set(gca,'xtick',1:4,'xticklabel',{'false','mixed','true','unclassified'})
        ylabel('# transients')
        xlabel('human classification')
        title('auto-classified as true')
        subplot(2,1,2)
        bar(Nfalse)
        set(gca,'xtick',1:4,'xticklabel',{'false','mixed','true','unclassified'})
        ylabel('# transients')
        xlabel('human classification')
        title('auto-classified as false')
    end
    
end





function [h,pValue,stat,cValue] = lbqtest2d(res,varargin)

% [h,pValue,stat,cValue] = lbqtest2d(res,varargin)
%
% 2-D extension of the Ljung-Box test
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input parsing

T             = sum(~isnan(res(:)));                                       % Effective sample size
nan_mask      = isnan(res);                                                % NaNs mask missing data
res(nan_mask) = 0;                                                         % Zero out missing data for correlation calculation

defaultLags  = [3,min(15,T-1)];                                            % Recommended in [1]
defaultAlpha = 0.05;                                                       % Default to rejection of the null at the 5% rate
defaultDof   = 0;                                                          % Store default lags

parseObj = inputParser;
parseObj.addRequired ('res'  ,               @(x) validateattributes(x, {'double'}, {'nonempty'}, '', 'residual series'));
parseObj.addParameter('Lags' , defaultLags , @(x) validateattributes(x, {'double'}, {'nonempty' 'integer'  'vector' 'positive' '<=' (T-1)}, '', 'number of lags'));
parseObj.addParameter('Alpha', defaultAlpha, @(x) validateattributes(x, {'double'}, {'nonempty' 'vector' '>' 0 '<' 1}, '', 'significance levels'));
parseObj.addParameter('DoF'  , defaultDof  , @(x) validateattributes(x, {'double'}, {'nonempty' 'vector' 'nonnegative' 'integer'}, '', 'degrees-of-freedom'));

try
    parseObj.parse(res,varargin{:});
catch ME
    throwAsCaller(ME)
end

res   = parseObj.Results.res;
lags  = parseObj.Results.Lags;
alpha = parseObj.Results.Alpha;

if any(strcmpi('DoF', parseObj.UsingDefaults))                             % If the user did not specify the degrees-of-freedom, then assign the number of lags as a default.
    dof = [];                                                              % Reset dof default
else
    dof = parseObj.Results.DoF;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the sample ACF out to the largest lag:

maxLag   = max(lags);                                                      % Get max lags
minLag   = min(lags);                                                      % Get min lags
wX       = size(res,1);                                                    % Get dimensions of residual image
wY       = size(res,2);                                                    % ---
[gX, gY] = meshgrid(-(wX-1):(wX-1), -(wY-1):(wY-1));                       % Create indexed grid of pixel locations
ACF      = xcorr2(nans2zeros(res));                                                    % Calculate autocorrelations
ACF      = ACF/ACF((gX==0)&(gY==0));                                       % Normalize the correlations to get coefficients
ACF      = ACF((sqrt(gX.^2+gY.^2)>minLag)&(sqrt(gX.^2+gY.^2)<maxLag));     % Select the correlations to test

num_var = xcorr2(double(~nan_mask));                                       % Get the number of variables used per correlation calculation
lat_vec = num_var((sqrt(gX.^2+gY.^2)>minLag)&(sqrt(gX.^2+gY.^2)<maxLag));  % Select the numbers corresponding to the used correlations

if isempty(dof)||(dof==0)
%     dof = numel(ACF);
    dof = sum(ACF~=0);
end
stat = T*(T+2)*nansum((ACF.^2)./lat_vec(:));                               % Compute Q-statistic

pValue = 1-chi2cdf(stat,dof);                                              % Compute p-values:

if nargout >= 4
    cValue = chi2inv(1-alpha,dof);                                         % Compute critical values, if requested:
else
    cValue = [];
end

h = (alpha >= pValue);                                                     % Perform the test:




function x = nans2zeros(x)
% one-liner to zero nans
x(isnan(x)) = 0;






