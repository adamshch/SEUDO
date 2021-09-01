function varargout = optimizeAutoClassification(se,varargin)
% results = seudoObject.optimizeAutoClassification(...)
%
% classify transients based on correlation between transient profile and source
% profile, then plot the ROC curve
%
% to meausre performance, ground truth is taken as the current classification 
%
%
% PARAMETERS
%   
%   threshValues - if scalar, interpreted as the number of thresholds to
%                       check in the range [-1, 1]
%                  if a vector, interpretted as the correlation values to check
%
%
% EXAMPLES
% 
%
%	se.optimizeAutoClassification(10)
%	se.optimizeAutoClassification(-.5:.05:.8)
%	results = se.optimizeAutoClassification;
%   
%
%
% 2018 - Adam Charles & Jeff Gauthier

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Input checking
%
% if ~iscell(se)
%     TMP   = se;
%     se    = cell(1);
%     se{1} = TMP;
%     clear TMP
% end
%
% if nargin > 1
%     n_pts = varargin{1};
% else
%     n_pts = 100;
% end


% parse inputs
p = inputParser;
p.addParameter('threshValues',20);
parse(p,varargin{:});


params = p.Results;
clear p

n_pts = params.threshValues;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get bounds to sweep over

min_corr = 0;
max_corr = 1;
for cc = 1:se.nCells
    min_corr = min(min_corr, ...
        min(se.tcDefault.transientInfo(cc).corrWithProfile));
    max_corr = max(max_corr, ...
        max(se.tcDefault.transientInfo(cc).corrWithProfile));
end

min_corr    = floor(10*min_corr)/10;                                       % Ensure that min correlation includes all transients
max_corr    = ceil(max_corr);                                              % Ensure that max correlation excludes all transients
if numel(n_pts) == 1
    if n_pts > 1
        thresh_vals = linspace(min_corr,max_corr,n_pts);                   % Set up the thresholds to test in the ROC plot
    else
        thresh_vals = min_corr:n_pts:max_corr;                             % Set up the thresholds to test in the ROC plot
        n_pts       = numel(thresh_vals);
    end
else
    thresh_vals = n_pts;
    n_pts       = numel(thresh_vals);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get ROC points
p_tp = zeros(n_pts,1);
p_fp = zeros(n_pts,1);

p_am = zeros(n_pts,1);
p_un = zeros(n_pts,1);

fprintf('Generating ROC curve')
for kk = 1:n_pts
    classHuman = [];
    classAuto  = [];
    CA = se.autoClassifyTransients('default', 'saveResults', ...
        false,'plotComparison',false, 'method', 'corr', ...
        'ignoreZeros',false, 'corrThresh', thresh_vals(kk)); % Calculate the threshold classification
    for cc = 1:se.nCells
        classHuman = cat(1,classHuman,...
            se.tcDefault.transientInfo(cc).classification); % Store the human classifications
        classAuto  = cat(1,classAuto,CA(cc).classification);           % Get the new automatic classifications
    end
    p_tp(kk) = sum((classHuman== seudo.valTrue)  & (classAuto==seudo.valTrue)) / sum(classHuman== seudo.valTrue);  	% Count the true positives
    p_fp(kk) = sum((classHuman== seudo.valFalse) & (classAuto==seudo.valTrue)) / sum(classHuman== seudo.valFalse); 	% Count the false positives
    p_am(kk) = sum((classHuman== seudo.valMix)   & (classAuto==seudo.valTrue)) / sum(classHuman== seudo.valMix);    % Count the mixed positives
    p_un(kk) = sum( isnan(classHuman)            & (classAuto==seudo.valTrue)) / sum(isnan(classHuman));            % Count the unclassified positives
    fprintf('.')
end
fprintf('\n')

h = figure;
set(h,'name',sprintf('%s, %d cells, %d transients, %d threshold values',se.name,se.nCells,length(classHuman),length(thresh_vals)))
subplot(2,3,[1,5])
axis image
xlim([0 1]);ylim([0 1])
set(gca,'xtick',0:.2:1,'ytick',0:.2:1,'tickdir','out')
grid on
hold on
box on
plot(p_fp,p_tp,'.-','color','k', 'linewidth',2,'markersize',20)
box off
axis square
set(gca,'TickDir','out','XTick',0:0.2:1,'YTick',0:0.2:1)
%xlabel('P(false positive)'); ylabel('P(true positive)')
xlabel('false'); ylabel('true'); title('fraction of transients labeled as true by auto-classifier')

% label data point closest to [0,1]
dists = (p_fp.^2) + ((1-p_tp).^2);
[~,closestPoint] = min(dists);
plotCol = [.5 0 1];
text(p_fp(closestPoint),p_tp(closestPoint),sprintf('\n  %0.2f',thresh_vals(closestPoint)),'fontsize',18,'color',plotCol)
plot(p_fp(closestPoint),p_tp(closestPoint),'.','color',plotCol,'markersize',30)

subplot(2,3,3)
plot(thresh_vals,p_tp,'color',seudo.pickColor('true'))
xlabel('correlation threshold')
ylabel('P(classsified as true)')
hold on
plot(thresh_vals,p_fp,'color',seudo.pickColor('false'))
box off
set(gca,'TickDir','out','YTick',0:0.5:1)

subplot(2,3,6)
plot(thresh_vals,p_am,'color',seudo.pickColor('mixed'))
hold on
plot(thresh_vals,p_un,'color',seudo.pickColor('unclassified'))
xlabel('correlation threshold')
ylabel('P(classsified as true)')
box off
set(gca,'TickDir','out','YTick',0:0.5:1)

set(gcf,'color',[1,1,1])


if nargout >=1
    extras.h     = h;
    extras.thresh_vals = thresh_vals;
    extras.p_fp  = p_fp;
    extras.p_tp  = p_tp;
    extras.p_am  = p_am;
    extras.p_un  = p_un;
    varargout{1} = extras;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%