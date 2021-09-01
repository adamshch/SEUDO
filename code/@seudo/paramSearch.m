function results = paramSearch(se,pickedCellTransients,factOrSets,varargin)
% results = seudoObject.paramSearch(pickedCellTransients,factOrSets,...)
%
% for a chosen set of transients, apply SEUDO using a variety of parameter combinations
%
%
% EXAMPLE
%
%   results = se.paramSearch([1 50 true; 2 16 false],'factorial','paramA',[1 2 3],'paramB',[0 0.5])
%
%
%
%
% INPUTS
%
%   pickedCellTransients - N x 3 matrix, specifies which transients will be used for the parameter search
%       each row is of the form [ cell number, transient number, true or false transient ]
%
%       e.g. seudoObject.paramSearch( [ 11 6 true ;  23 2 false ;  78 5 false ] ,...)
%
%       (specifying logical values improves readability, though these values are converted to double in the matrix)
%
%
%   factOrSets - can be 'factorial', 'sets', or a cell array of structs
%
%   subsequent arguments specify the parameters to search in one of three ways:
%   
%   1) factorial design, with vectors or cell arrays giving the values for each parameter
%
%       for example, this input:
%
%           paramSearch(...,'factorial','paramA',[1 2 3],'paramB',{'up','down'})
%
%       results in executing these SEUDO commands:
%
%           estimateTimeCoursesWithSEUDO('paramA',1,'paramB','up')
%           estimateTimeCoursesWithSEUDO('paramA',1,'paramB','down')
%           estimateTimeCoursesWithSEUDO('paramA',2,'paramB','up')
%           estimateTimeCoursesWithSEUDO('paramA',2,'paramB','down')
%           estimateTimeCoursesWithSEUDO('paramA',3,'paramB','up')
%           estimateTimeCoursesWithSEUDO('paramA',3,'paramB','down')
%
%
%   2) sets, with vectors or cell arrays of possible parameter values, each of length 1 or N
%
%       for example, this input:
%
%           paramSearch(...,'sets','paramA',[1 2 3],'paramB',{'up','down','off'},'paramC',0)
%
%       results in executing these SEUDO commands:
%
%           estimateTimeCoursesWithSEUDO( 'paramA',1, 'paramB','up'  , 'paramC',0 )
%           estimateTimeCoursesWithSEUDO( 'paramA',2, 'paramB','down', 'paramC',0 )
%           estimateTimeCoursesWithSEUDO( 'paramA',3, 'paramB','off' , 'paramC',0 )
%
%
%
%   3) a cell array of structs, each fully specifying the parameters to use
%
%       for example, this input:
%
%           paramSearch(...,{ struct('paramA',1,'paramB','off') , struct('paramA',2,'paramC',0) })
%
%       results in executing these SEUDO commands:
%
%           estimateTimeCoursesWithSEUDO('paramA',1,'paramB','off')
%           estimateTimeCoursesWithSEUDO('paramA',2,'paramC',0)
%
%   
%
%
%
%   for numeric values, a vector is treated as a list of possible parameter values
%   to specify a vector as a parameter value, provide it in cell string, e.g.
%
%       estimateTimeCoursesWithSEUDO('paramA', [1 2 3] )     is interpretted as THREE possible SCALAR values for paramA
%       estimateTimeCoursesWithSEUDO('paramA',{[1 2 3]})     is interpretted as ONE possible VECTOR value for paramA
%
%
%
%
% OUTPUTS
%
%   results - struct containing results of parameter search, should be used as an argument to seudoReviewParamSearch
%
%
% 2018 Jeff Gauthier



% VERIFY INPUTS
% (more verification below)

if size(pickedCellTransients,2) ~= 3
    error('arugment ''pickedCellTransients'' must have 3 columns')
end





% MAKE LIST OF PARAMETER SETS



% if a cell array was provided, assume it is a cell array of structs, each with the values for one parametet set
if iscell(factOrSets)
    
    paramSets = factOrSets;
    
    % ensure all are structs
    if ~all(cellfun(@isstruct,paramSets))
        error('second argument must be ''factorial'', ''sets'', or a cell array of structs')
    end
    
else
    
        
    
    
    
    % ensure even
    if mod(length(varargin),2) ~= 0
        error('Second arugment must be either ''factorial'', ''sets'', or a cell array of structs. Try ''help paramSearch''.')
    end
    
    % get counts
    nParams         = length(varargin)/2;
    paramCounts     = cellfun(@length,varargin(2:2:end));
    singletonParams = paramCounts == 1;
    
    paramNames             = varargin(1:2:end);
    singletonParamNames    = paramNames(singletonParams);
    nonSingletonParamNames = paramNames(paramCounts~=1);
    
    
    % ensure no duplicate field names
    if length(paramNames) ~= length(unique(paramNames))
        [uniqueNames,~,whichSpots] = unique(paramNames);
        n = hist(whichSpots,1:max(whichSpots));
        repName = uniqueNames{find(n>1,1)};
        error('The parameter ''%s'' was repeated. All parameters must be specified only once.',repName)
    end
    
    
    
    % GENERATE PARAMETER SETS
    
    % create cell array of structs, each containing the parameter specs for one run
    
    
    switch factOrSets
        case 'factorial'
            
            
            % design full factorial analysis
            paramCodex = fullfact(paramCounts);
            
            % populate parameter sets
            paramSets = cell(size(paramCodex,1),1);
            
            
            for pp = 1:nParams
                
                paramName   = varargin{2*pp-1};
                paramValues = varargin{2*pp};
                
                % if numeric, sort parameter values 
                if isnumeric(paramValues)
                    paramValues = sort(paramValues);
                end
                
                % put values into the field for each parameter set
                for ss = 1:size(paramCodex,1)
                    paramSets{ss}.(paramName) = getElement(paramValues,paramCodex(ss,pp));
                end
            end
            
            
        case 'sets'
            
            
            % ensure all non-singleton params all have the same length
            if length(unique(paramCounts(~singletonParams))) ~= 1
                error(['In ''sets'' specification, each parameter must have either 1 argument or N arguments. ' ...
                    'The provided parameters had multiple lengths: %s'],num2str(paramCounts))
            end
            
            
            nSets = unique(paramCounts(~singletonParams));
            
            paramSets = cell(nSets,1);
            for ss = 1:nSets
                for pp = 1:nParams
                    
                    % choose which element to get from the list of possible parameter values
                    if singletonParams(pp)
                        indToGet = 1;
                    else
                        indToGet = ss;
                    end
                    
                    % get parameter values
                    paramValues = varargin{2*pp};
                    
                    % put into this parameter set
                    paramSets{ss}.(varargin{2*pp-1}) = getElement(paramValues,indToGet);
                end
            end
            
            
        otherwise
            error('second argument must be ''factorial'', ''sets'', or a cell array of structs')
            
    end
    
end




% generate paramCodex, a matrix that indicates which value of each parameter is used in each parameter set
%   rows are parameter sets, columns are parameters
%   for a parameter with N possible values:
%       an integer k in [1,N] indicates the k-th possible value
%       0 means the field does not occur so the default value is used
%
if ~exist('paramCodex','var')
    
    paramsPerSet = cellfun(@fieldnames,paramSets,'uniformoutput',0);
    
    paramNames = unique(cat(1,paramsPerSet{:}));
    
    paramCodex = zeros(length(paramSets),length(paramNames));
    
    % identify possible values for each parameter
    paramValsOfEach   = cell(length(paramNames),1);
    possibleParamVals = cell(length(paramNames),1);
    for pp = 1:length(paramNames)
        setsWithThisParam     = cellfun(@(x)isfield(x,paramNames{pp}),paramSets);
        paramValsOfEach{pp}   = cellfun(@(x)x.(paramNames{pp}),paramSets(setsWithThisParam),'uniformoutput',0);
        possibleParamVals{pp} = getUniqueValues(paramValsOfEach{pp});
    end
    
    
    for ss = 1:length(paramSets)
        
        for pp = 1:length(paramNames)
            
            if isfield(paramSets{ss},paramNames{pp})
                %paramCodex(ss,pp) = find(ismember(paramValsOfEach{pp},{paramSets{ss}.(paramNames{pp})}));
                paramCodex(ss,pp) = find(cellfun(@(x)isequaln(x,paramSets{ss}.(paramNames{pp})),possibleParamVals{pp}));
            end
        end
    end
    
    % identify singleton, non-singleton params
    for pp=1:length(paramNames)
        paramCounts(pp) = length(unique(paramCodex(:,pp)));
    end
    singletonParamNames    = paramNames(paramCounts == 1);
    nonSingletonParamNames = paramNames(paramCounts~=1);
end








% PERFORM ANALYSIS ON EACH PARAMETER SET




% put each transient into its own seudo object
minis = makeMiniSeudos(se,pickedCellTransients(:,1:2),'allCells',true);



% ensure all parameter sets are valid before executing any of them
for mm = 1:length(minis)
    for ss = 1:length(paramSets)
        res = minis{mm}.estimateTimeCoursesWithSEUDO('validateParams',true,paramSets{ss});
        if ~isempty(res)
            fprintf('\n************************************\n     PROBLEM IN PARAMETER SET %d    \n************************************\n\n',ss)
            rethrow(res)
        end
    end
end





% show progress
startTime = now;
fprintf('dataset %s:\n   will perform SEUDO for %d sets of parameters on %d transients (%d total)...\n',...
    se.name,length(paramSets),length(minis),length(paramSets)*length(minis))


% pick on which seudo calculations to show update
nSeudo      = length(minis) * length(paramSets);
initialSets = [1 3 10 30 100 300 1000];
initialSets = initialSets(initialSets<nSeudo*.2);
prcSets     = round(linspace(0,nSeudo,21));
if ~isempty(initialSets)
    prcSets     = prcSets(prcSets>initialSets(end));
    updateSets  = [initialSets prcSets];
else
    updateSets  = prcSets;
end
    



tcResult = struct;

% cycle through all transients (one per seudo object)
for mm = 1:length(minis)
    
    % cycle through parameter sets
    for ss = 1:length(paramSets)
        
        % perform SEUDO
        minis{mm}.estimateTimeCoursesWithSEUDO('verbose',0,paramSets{ss},'whichCells',1);
        
        % provide progress update
        nn = (mm-1)*length(paramSets) + ss;
        if ismember(nn,updateSets)
            finishDur = (nSeudo-nn)/nn*(now-startTime);
            fprintf('      completed %d of %d (%0.1f%%) in %0.1f min, estimated finish in %0.1f min (%s)\n',...
                nn, nSeudo, 100*nn/nSeudo,(now-startTime)*24*60,finishDur*24*60,datestr(now+finishDur))
        end
        
    end
    
    % store result
    tcResult(mm).tcDefault = minis{mm}.tcDefault;
    tcResult(mm).tcSeudo = minis{mm}.tcSeudo;
    
    
    
end

fprintf('finished in %0.1f min\n',(now-startTime)*24*60)






% PACKAGE RESULTS

results                        = struct;
results.pickedCellTransients   = pickedCellTransients;
results.tc                     = tcResult;
results.paramSets              = paramSets;
results.paramCodex             = paramCodex;
results.paramNames             = paramNames;
results.singletonParamNames    = singletonParamNames;
results.nonSingletonParamNames = nonSingletonParamNames;
results.argspecs               = varargin;
results.name                   = se.name;
results.startTime              = startTime;




function uv = getUniqueValues(input)
%identify unique elements in a cell array (matlab can't do this in general)

uv = input(1);

for cc = 2:length(input)
    if ~any(cellfun(@(x)isequaln(input{cc},x),uv))
        uv{end+1} = input{cc};
    end
end



function element = getElement(array,N)
% get Nth element in an array (either vector or cell)

if iscell(array)
    element = array{N};
else
    element = array(N);
end

