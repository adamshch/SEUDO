function results = bayesianParameterOptimization(se,transSpec,varargin)
% seudoObject.bayesianParameterOptimization(transSpec)
%
% select optimal parameters using Bayesian optimization
%
%   transSpec - T x 2 matrix specifying a set of transients
%               each row is [ < cell number >  < transient number > ]
%
%
% 2018 Adam Charles

p = inputParser;
p.addParameter('p',[]);                                                    % probability of only the neural profile
p.addParameter('sigma2',[]);                                               % noise variants
p.addParameter('lambdaBlob',[]);                                           % cost per Gaussian blob
p.addParameter('blobRadius',[]);                                           % Size of the Gaussian blobs
parse(p,varargin{:});


defaults.p          = [0,1];
defaults.sigma2     = [0,200];
defaults.lambdaBlob = [1e-4,2];
defaults.blobRadius = [0.5,4.5];

paramNames  = fieldnames(p.Results);
vars_to_opt = [];
eval_str    = 'fun = @(x) computeError(se.makeMiniSeudos(transSpec),transSpec(:,3)';


for kk = 1:numel(paramNames)
    eval(['TMP = p.Results.',paramNames{kk},';']);
    if any(isnan(TMP))||isempty(TMP)                                       % If empty or a nan, optimize this parameter
        if isfield(defaults,paramNames{kk})                                % Check if a default range for this parameter
            eval(['TMP2 = defaults.',paramNames{kk},';']);                 % Check if a default range for this parameter
            vars_to_opt = [vars_to_opt,...
                                optimizableVariable(paramNames{kk},TMP2)]; % create an optimization variable over the default range
        else
            vars_to_opt = [vars_to_opt,...
                          optimizableVariable(paramNames{kk},[1e-5,100])]; % create an optimization variable over the default default range
        end
        eval_str = [eval_str,',''',paramNames{kk},''',x.',paramNames{kk}]; % Add the parameter appropriately to the anonymous function evaluation string
    elseif numel(TMP)==2                                                   % If a range, optimize this parameter over that range
        vars_to_opt(kk) = [vars_to_opt,...
                 optimizableVariable(paramNames{kk},[min(TMP),max(TMP)])]; % create an optimization variable over the given range
        eval_str = [eval_str,',''',paramNames{kk},''',x.',paramNames{kk}]; % Add the parameter appropriately to the anonymous function evaluation string
    elseif numel(TMP)==1
        eval_str = [eval_str,',''',paramNames{kk},''',x.',...
                                                       sprintf('%f',TMP)]; % If a single value, pass that on to SEUDO
    else                                                                   % Use default parameters
    end
end

eval_str = [eval_str,');'];                                                % Cap off the evaluation string
eval(eval_str);                                                            % Evaluate the string to create the function to optimize

% perform optimization
results = bayesopt(fun,vars_to_opt,'Verbose',1,'AcquisitionFunctionName','expected-improvement-plus','IsObjectiveDeterministic',true);


function netError = computeError(minis,transClass,varargin)
% compute error


% do SEUDO for each
nTransients = length(minis);
for ss = 1:nTransients
    minis{ss}.estimateTimeCoursesWithSEUDO(varargin{:},'verbose',0);
end


% compute what fraction of each transient was still present after running SEUDO
keptFractions = nan(nTransients,1);
for ss = 1:nTransients
    
    % original time course
    origTrace = minis{ss}.tcDefault.tc;
    
    % seudo time course
    seudoTrace = minis{ss}.tcSeudo(end).tc;
    
    % fraction still present
    keptFractions(ss) = sum(seudoTrace)/sum(origTrace);
end


% use provided classification to evaluate performance
meanTrue = nanmean(keptFractions(transClass==1));
meanFalse = nanmean(keptFractions(transClass==0));


% combine in a simple way 
netError = meanFalse - meanTrue;



