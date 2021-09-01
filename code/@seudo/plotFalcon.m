function varargout = plotFalcon(se,varargin)
% Generate Frequencey and Level of Contamination (FaLCon) plot
%
% seudoObject.getFalconPlotPoints(...)
%
% for full usage and options, see help for function getFalconPlotPoints
% (inputs to plotFalcon are directly passed to getFalconPlotPoints)
%
%
% EXAMPLES:
%
%   seudoObject.getFalconPlotPoints 
%   seudoObject.getFalconPlotPoints([],[],'alpha',0.03,'radius',[2 7])
%
%
%
%
% 2018 Jeff Gauthier


[x,y] = se.getFalconPlotPoints(varargin{:});

countPerCell = [];
for cc = 1:se.nCells
    C             = se.tcDefault.transientInfo(cc).classification;
    countPerCell  = cat(1, countPerCell, [sum(C==1) sum(C==0) sum(C==-1), sum(isnan(C))]);
end
h_falcon = plotFalconInnerFcn(x,y, 'counts', sum(countPerCell,2), 'numTrue', countPerCell(:,1), 'numClass', sum(countPerCell(:,1:3),2));



if nargout > 0
    varargout{1} = h_falcon;
end


end




function h = plotFalconInnerFcn(m_cfrac, m_rfrac, varargin)

% varargout = plotFaLCon(m_cfrac,m_rfrac, varargin)
%
% Generate FaLCon plot
%  - m_cfrac - x-axis of points in the FaLCon plot
%  - m_rfrac - y-axis of points in the FaLCon plot
%  - Extra parameters (call using ... 'paramName', paramValue...)
%   - counts       (default = [])    - The counts of how many transients in
%                                      a point 
%   - fracs        (default = [])    - The fraction of true transients in
%                                      the point 
%   - figNum       (default = [])    - The figure to plot to
%   - fracBounds   (default = [])    - Bounds of fraction true to restrict
%                                      plotting to 
%   - dispProfNums (default = false) - Toggle whether to display cell
%                                      numbers with the points 
%   - nBins        (default = 20)    - number of bins for the histogram 
%
% Output:
%  - h - handle of the figure, or array of figure handles if multiple plots
%        were generated
%
% 2018 Adam Charles & Jeff Gauthier

p = inputParser;
p.addParameter('counts'      , []);                                        % The counts of how many transients in a point
p.addParameter('fracs'       , []);                                        % The fraction of true transients in the point

p.addParameter('numTrue'     , []);                                        % The number of true transients in the point
p.addParameter('numClass'    , []);                                        % The number of transients classified for each cell


p.addParameter('figNum'      , []);                                        % The figure to plot to
p.addParameter('fracBounds'  , []);                                        % Bounds of fraction true to restrict plotting to
p.addParameter('dispProfNums', false);                                     % Toggle whether to display cell numbers with the points
p.addParameter('nBins'       , 20);                                        % number of bins for the histogram 
p.addParameter('zeroJitter'  , true);                                      % toggle whether or not to jitter points at zero
parse(p,varargin{:});

params = p.Results;
clear p

if params.zeroJitter
    IDX_allzeros          = (m_cfrac==0)&(m_rfrac==0);
    m_cfrac(IDX_allzeros) = 0.01*rand(sum(IDX_allzeros),1);
    m_rfrac(IDX_allzeros) = 0.01*rand(sum(IDX_allzeros),1);
end

if (~isempty(params.numClass))&&(~isempty(params.numTrue))
    params.fracs = params.numTrue./params.numClass;
end

if isempty(params.counts)
    IDX = 1:numel(m_cfrac);
    h   = do_plotting(m_cfrac, m_rfrac, [], [], IDX, params.nBins,...
                                      params.dispProfNums, params.figNum);
else
    if any(isnan(params.fracBounds))
        IDX = find(params.counts==0);
        h   = do_plotting(m_cfrac, m_rfrac, params.counts, params.fracs, ...
                     IDX, params.nBins, params.dispProfNums, params.figNum);
    elseif isnumeric(params.fracBounds)&&(numel(params.fracBounds)==2)
        IDX = find((params.fracs>min(params.fracBounds))&(params.fracs<max(params.fracBounds)));
        h   = do_plotting(m_cfrac, m_rfrac, params.counts, params.fracs, ...
                     IDX, params.nBins, params.dispProfNums, params.figNum);
    elseif isnumeric(params.fracBounds)&&(numel(params.fracBounds)==1)
        params.fracBounds = [params.fracBounds round(params.fracBounds)];
        IDX = find((params.fracs>min(params.fracBounds))&(params.fracs<max(params.fracBounds)));
        h   = do_plotting(m_cfrac, m_rfrac, params.counts, params.fracs, IDX, ...
                         params.nBins, params.dispProfNums, params.figNum);
    elseif isnumeric(params.fracBounds)&&(size(params.fracBounds,2)==2)&&(size(params.fracBounds,1)>1)
        h = zeros(size(params.fracBounds,1),1);                            % Initialize the graphical object array
        for ll = 1:size(params.fracBounds,1)
            h(ll) = plotFalconInnerFcn(m_cfrac, m_rfrac, 'counts', params.counts, ...
                   'fracs', params.fracs, 'fracBounds', params.fracBounds(ll,:), ...
                   'nBins', params.nBins, 'dispProfNums', params.dispProfNums, 'figNum', params.figNum); 
        end
    elseif isnumeric(params.fracBounds)&&(size(params.fracBounds,1)==2)&&(size(params.fracBounds,2)>1)
        h = zeros(size(params.fracBounds,1),1);                            % Initialize the graphical object array
        for ll = 1:size(params.fracBounds,2)
            h(ll) = plotFalconInnerFcn(m_cfrac,m_rfrac, 'counts', params.counts, ...
                   'fracs', params.fracs, 'fracBounds', params.fracBounds(:,ll), ...
                   'nBins', params.nBins, 'dispProfNums', params.dispProfNums, 'figNum', params.figNum);
        end
    else
        IDX = 1:numel(m_cfrac);
        h = do_plotting(m_cfrac ,m_rfrac, params.counts, params.fracs, IDX, ...
                         params.nBins, params.dispProfNums, params.figNum);
    end
end

end

function h = do_plotting(m_cfrac,m_rfrac,countPerCell,TMPfrac,IDX,n_bins,dispProfNums,figno)

if ~isempty(figno)
    h = figure(figno);
else
    h = figure();
end

base_diam = 10;

ix      = false(numel(m_cfrac),1);
ix(IDX) = true;
IDX     = ix;
clear ix

countPerCell(countPerCell<10) = nan;

subplot(3,3,[1,2]); 
histogram(m_cfrac(IDX),linspace(0,1,n_bins));
set(gca, 'XLim',[0,1], 'XTick', [0,1])
ylabel('# sources')
top_axs = gca;
box off

if any(isnan(TMPfrac))||any(isnan(countPerCell))||any(countPerCell==0)
    
    IDX_nan   = (vec(IDX)&vec(isnan(TMPfrac(IDX)))|(countPerCell(IDX)==0));
    IDX_nonan = vec(IDX)&vec((~isnan(TMPfrac(IDX)))&(countPerCell(IDX)>0));
    
    subplot(3,3,[4,8]), cla
    subplot(3,3,[4,8]), hold on
    
    if isempty(countPerCell)
%         subplot(3,3,[4,8]); scatter(m_cfrac(IDX_nan),m_rfrac(IDX_nan),base_diam,'k','filled');
        subplot(3,3,[4,8]); scatter(m_cfrac(IDX_nan),m_rfrac(IDX_nan),base_diam,'k');
    else
        subplot(3,3,[4,8]); scatter(m_cfrac(IDX_nan),m_rfrac(IDX_nan),base_diam+sum(countPerCell(IDX_nan,:),2),'k');
%         subplot(3,3,[4,8]); scatter(m_cfrac(IDX_nan),m_rfrac(IDX_nan),base_diam+sum(countPerCell(IDX_nan,:),2),'k','filled');
    end
    
    if (~isempty(countPerCell))&&isempty(TMPfrac)
        subplot(3,3,[4,8]); scatter(m_cfrac(IDX_nonan),m_rfrac(IDX_nonan),base_diam+countPerCell(IDX_nonan,:),[],'filled');
    elseif isempty(countPerCell)&&(~isempty(TMPfrac))
        subplot(3,3,[4,8]); scatter(m_cfrac(IDX_nonan),m_rfrac(IDX_nonan),base_diam,TMPfrac(IDX_nonan),'filled');
    else
        subplot(3,3,[4,8]); scatter(m_cfrac(IDX_nonan),m_rfrac(IDX_nonan),base_diam+sum(countPerCell(IDX_nonan,:),2),TMPfrac(IDX_nonan),'filled');
    end
    

    subplot(3,3,[4,8]), hold off
   
else
    if isempty(countPerCell)&&isempty(TMPfrac)
        subplot(3,3,[4,8]); scatter(m_cfrac(IDX),m_rfrac(IDX),base_diam,'filled')
    elseif (~isempty(countPerCell))&&isempty(TMPfrac)
        subplot(3,3,[4,8]); scatter(m_cfrac(IDX),m_rfrac(IDX),base_diam+countPerCell(IDX),'filled');
    elseif isempty(countPerCell)&&(~isempty(TMPfrac))
        subplot(3,3,[4,8]); scatter(m_cfrac(IDX),m_rfrac(IDX),base_diam,TMPfrac(IDX),'filled');
    else
        subplot(3,3,[4,8]); scatter(m_cfrac(IDX),m_rfrac(IDX),base_diam+sum(countPerCell(IDX,:),2),TMPfrac(IDX),'filled');
    end
end

if dispProfNums
    subplot(3,3,[4,8]), hold on
    for tt=1:length(m_cfrac)
        text(double(m_cfrac(tt)),double(m_rfrac(tt)),num2str(tt))          % Display the number for the tt^th cell on the scatter plot
    end
    subplot(3,3,[4,8]), hold off
end
caxis([0 1])
cmap = [linspace(1,0,100)' zeros(100,1) linspace(0,1,100)'];
colormap(cmap)
mid_axs = gca;
set(gca, 'XLim', [0,1], 'YLim', [0,1])
axis square
xlabel('Contaminated rate')
ylabel('Contamination severity')

subplot(3,3,[6,9]); histogram(m_rfrac(IDX),linspace(0,1,n_bins));
set(gca, 'XLim',[0,1], 'XTick', [0,1])
ylabel('# sources')
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

if (~isempty(countPerCell))||(~isempty(TMPfrac))
    subplot(3,3,3), cla
    subplot(3,3,3), hold on
    
    if ~isempty(countPerCell)
        leg_diams = base_diam+linspace(min(countPerCell),max(countPerCell),8);
        scatter(vec(1:8),0.5*ones(8,1),vec(leg_diams),'k','filled');
        text(1,1,num2str(double(leg_diams(1))))
        text(8,1,num2str(double(leg_diams(end))))
    end
    if ~isempty(TMPfrac)
        cmap = [linspace(1,0,100)' zeros(100,1) linspace(0,1,100)'];
        colormap(cmap)
        cb = colorbar;
        cb.Location = 'northoutside';
        cb.Label.String = '% true';
    end
    axis off
    subplot(3,3,3), hold off
end


end


function Y = vec(X)
Y = X(:);
end