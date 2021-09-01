function varargout = plotTransientFaLCon(se, varargin)

% varargout = plotTransientFaLCon(m_cfrac,m_rfrac, varargin)
%
% Generate FaLCon plot
%  - m_cfrac - x-axis of points in the FaLCon plot
%  - m_rfrac - y-axis of points in the FaLCon plot
%  - Extra parameters (call using ... 'paramName', paramValue...)
%   - classType    - Classification of the transients (default to using classification stored in seudo object) 
%   - sepScatter   - toggle whether or not to put different classes of data
%                    on different plots (default = false)
%   - resSelect    - toggle what type of ratio to show on the y-axis
%                    (default = 'ratio') 
%   - figNum       - The figure to plot to (default = [])
%   - dispProfNums - Toggle whether to display cell numbers with the points
%                    (default = false) 
%   - nBins'       - number of bins for the histogram (default = 20)
%
% Output:
%  - h - handle of the figure, or array of figure handles if multiple plots
%        were generated
%
% 2018 Adam Charles & Jeff Gauthier

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input parsing

p = inputParser;
p.addParameter('figNum'      , []);                                        % The figure to plot to
p.addParameter('classType'   , []);                                        % Array of how data was classified
p.addParameter('dispProfNums', false);                                     % Toggle whether to display cell numbers with the points
p.addParameter('nBins'       , 20);                                        % number of bins for the histogram 
p.addParameter('zeroJitter'  , true);                                      % toggle whether or not to jitter points at zero
p.addParameter('sepScatter'  , false);                                     % toggle whether or not to put different classes of data on different plots 
p.addParameter('resSelect'   , 'ratio');                                   % toggle what type of ratio to show on the y-axis
parse(p,varargin{:});
params = p.Results;
clear p

if ~iscell(se)
    TMP   = se;
    se    = cell(1);
    se{1} = TMP;
    clear TMP
end

if ~(isequal(params.resSelect, 'ratio')||isequal(params.resSelect, 'norms')||isequal(params.resSelect, 'noCell'))
    warning('Unknown y-axis selection; reverting to plotting the residual ratio.')
    params.resSelect = 'ratio';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get transient Info

resRatios = [];                                                            % Initialize the residual ratio values (y-axis)
corrVals  = [];                                                            % Initialize the correlation values (x-axis)
classType = [];                                                            % Initialize the classification values (color)

for ll = 1:numel(se)
    for cc = 1:size(se{ll}.tcDefault.tc,2)
        tc_corrs  = se{ll}.tcDefault.transientInfo(cc).corrWithProfile;    % Extract correlation information from SEUDO cell
        r         = calcResRatios(se{ll}.profiles, ...
                                  se{ll}.tcDefault.transientInfo(cc), cc); % Extract the residual information from the seudo cell
        if isequal(params.resSelect, 'ratio')
            resRatios = cat(1,resRatios,vec(r.resRatios));                 % Concatenate the residual ratio values to the full vector
        elseif isequal(params.resSelect, 'norms')
            resRatios = cat(1,resRatios,vec(r.resNorms));                  % Concatenate the residual ratio values to the full vector
        elseif isequal(params.resSelect, 'noCell')
            resRatios = cat(1,resRatios,vec(r.resNoCell));                 % Concatenate the residual ratio values to the full vector
        end
        corrVals  = cat(1,corrVals,vec(tc_corrs));                         % Concatenate the correlation values to the full vector
        if isempty(params.classType)&&isfield(se{ll}.tcDefault.transientInfo(cc),'classification')
            classType = cat(1,classType,vec(se{ll}.tcDefault.transientInfo(cc).classification));
        end
    end
end

if max(resRatios(:))>1
    resRatios = resRatios/max(resRatios(:));                               % Normalize the data to [0,1] if a non-normalized metric was chosen
end

if isempty(params.classType)
    params.classType = classType;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Take care of plotting

if (params.sepScatter)&&(numel(params.classType)==numel(corrVals))
    
    trueOnes  = classType==seudo.valTrue;
    falseOnes = classType==seudo.valFalse;
    uncOnes   = isnan(classType);
    mixOnes   = ~isnan(classType) & ~ismember(classType,[seudo.valTrue seudo.valFalse]);
    
    h(1) = do_plotting(corrVals(trueOnes), resRatios(trueOnes), params.nBins, ...       % true
                    params.classType(trueOnes), params.dispProfNums, params.figNum);
                
    h(2) = do_plotting(corrVals(mixOnes), resRatios(mixOnes), params.nBins, ...         % mixed
                    params.classType(mixOnes), params.dispProfNums, params.figNum);
                
    h(3) = do_plotting(corrVals(falseOnes), resRatios(falseOnes), params.nBins, ...     % false
                    params.classType(falseOnes), params.dispProfNums, params.figNum);
                
    h(4) = do_plotting(corrVals(uncOnes), resRatios(uncOnes), params.nBins, ...         % unclassified
                    params.classType(uncOnes), params.dispProfNums, params.figNum);
else
    h    = do_plotting(corrVals, resRatios, params.nBins, ...
                    params.classType, params.dispProfNums, params.figNum);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output plotting

if nargout > 0
    varargout{1} = h;
end

end

function h = do_plotting(corrVals, resRatios, n_bins, classType, dispProfNums,figno)

if ~isempty(figno)
    h = figure(figno);
else
    h = figure();
end

subplot(3,3,[1,2]); histogram(corrVals,linspace(0,1,n_bins));
set(gca, 'XLim',[0,1])
top_axs = gca;
box off

if (~isempty(classType))&&(numel(classType)==numel(corrVals))
    
    trueOnes  = classType==seudo.valTrue;
    falseOnes = classType==seudo.valFalse;
    uncOnes   = isnan(classType);
    mixOnes   = ~isnan(classType) & ~ismember(classType,[seudo.valTrue seudo.valFalse]);
    
    subplot(3,3,[4,8]); scatter(corrVals(mixOnes), resRatios(mixOnes),'g','filled');
    subplot(3,3,[4,8]), hold on
    subplot(3,3,[4,8]); scatter(corrVals(trueOnes), resRatios(trueOnes),'b','filled');
    subplot(3,3,[4,8]); scatter(corrVals(falseOnes),resRatios(falseOnes),'r','filled');
    subplot(3,3,[4,8]); scatter(corrVals(uncOnes),resRatios(uncOnes),'k','filled');
    subplot(3,3,[4,8]), hold off
else
    subplot(3,3,[4,8]); scatter(corrVals,resRatios,'k','filled');
end

if dispProfNums
    subplot(3,3,[4,8]), hold on
    for tt=1:length(corrVals)
        text(corrVals(tt),resRatios(tt),num2str(tt))
    end
    subplot(3,3,[4,8]), hold off
end
caxis([0 1])
cmap = [linspace(1,0,100)' zeros(100,1) linspace(0,1,100)'];
colormap(cmap)
mid_axs = gca;
set(gca, 'XLim', [0,1], 'YLim', [0,1])
axis square
xlabel('Correlation value')
ylabel('Residual ratio')

subplot(3,3,[6,9]); histogram(resRatios,linspace(0,1,n_bins));
set(gca, 'XLim',[0,1])
view(90, -90)
box off

% Find current position [x,y,width,height]
pos1 = get(top_axs, 'Position');
pos2 = get(mid_axs, 'Position'); 

% Set width of second axes equal to first
pos1(1) = pos2(1)+0.1;
pos1(3) = pos2(3)-0.2;
set(top_axs,'Position',pos1);
set(gcf,'color',[1,1,1])

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%