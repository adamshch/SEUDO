function varargout = plotLBQDetectROC(se,varargin)

% [h] = plotCorrDetectROC(se,[n_pts])
% 
% plot the ROC curve using the correlation-based transient classification
% (potential x-axis for the FaLCon plots).
% 
% 2018 - Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input checking
if ~iscell(se)
    TMP   = se;
    se    = cell(1);
    se{1} = TMP;
    clear TMP
end

if nargin > 1
    n_pts = varargin{1};
else
    n_pts = 100;
end

if nargin > 2
    alpha_vals = varargin{2};
else
    alpha_vals = logspace(-9,-1,5);
end

num_alpha = numel(alpha_vals);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get bounds to sweep over

if n_pts > 1
    thresh_vals = linspace(0,1,n_pts);                                     % Set up the thresholds to test in the ROC plot
else
    thresh_vals = 0:n_pts:1 + eps;                                         % Set up the thresholds to test in the ROC plot
    n_pts       = numel(thresh_vals);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get ROC points
p_tp = zeros(n_pts,num_alpha);
p_fp = zeros(n_pts,num_alpha);

p_am = zeros(n_pts,num_alpha);
p_un = zeros(n_pts,num_alpha);

fprintf('Generating ROC curve...\n')
for aa = 1:num_alpha
    for kk = 1:n_pts
        fprintf('calculating for alpha = %f, thresh = %f', alpha_vals(aa),thresh_vals(kk))
        classHuman = [];
        classAuto  = [];
        for ll = 1:numel(se)
            CA = se{ll}.autoClassifyTransients('default','saveResults',...
                          false,'plotComparison',false,'method', 'lbq',...
                             'ignoreZeros',true,'radius',[4 6],'alpha',...
                         alpha_vals(aa),'tauBounds',thresh_vals(kk)*[1,1]);% Calculate the threshold classification
            for cc = 1:se{ll}.nCells
                classHuman = cat(1,classHuman,...
                           se{ll}.tcDefault.transientInfo(cc).classification); % Store the human classifications
                classAuto  = cat(1,classAuto,CA(cc).classification);           % Get the new automatic classifications
            end
            fprintf('.')
        end
        p_tp(kk,aa) = sum((classHuman== seudo.valTrue ) & (classAuto==seudo.valTrue)) / sum(classHuman== seudo.valTrue);     % Count the true positives
        p_fp(kk,aa) = sum((classHuman== seudo.valFalse) & (classAuto==seudo.valTrue)) / sum(classHuman== seudo.valFalse);    % Count the false positives
        p_am(kk,aa) = sum((classHuman== seudo.valMix  ) & (classAuto==seudo.valTrue)) / sum(classHuman== seudo.valMix);      % Count the mixed positives
        p_un(kk,aa) = sum( isnan(classHuman)            & (classAuto==seudo.valTrue)) / sum(isnan(classHuman));              % Count the unclassified positives
        fprintf('\n')
    end
end

h = figure;
subplot(2,4,[1,6]), plot(p_fp(:,1),p_tp(:,1),'g')
subplot(2,4,[1,6]), hold on
for ll = 2:size(p_fp,2)
    subplot(2,4,[1,6]), plot(p_fp(:,ll),p_tp(:,ll),'g')
end
xlabel('P(false positive)')
ylabel('P(true positive)')
box off
axis square
set(gca,'TickDir','out','XTick',0:0.2:1,'YTick',0:0.2:1)
subplot(2,4,[1,6]), hold off

subplot(2,4,3), plot(thresh_vals,p_tp,'b')
xlabel('Correlation threshold')
ylabel('P(true positive)')
box off
set(gca,'TickDir','out','YTick',0:0.5:1)

subplot(2,4,4), plot(thresh_vals,p_fp,'r')
xlabel('Correlation threshold')
ylabel('P(false positive)')
box off
set(gca,'TickDir','out','YTick',0:0.5:1)

subplot(2,4,7), plot(thresh_vals,p_am,'g')
xlabel('Correlation threshold')
ylabel('P(mixed positive)')
box off
set(gca,'TickDir','out','YTick',0:0.5:1)

subplot(2,4,8), plot(thresh_vals,p_un,'color', [0.5,0.5,0.5])
xlabel('Correlation threshold')
ylabel('P(uclassified positive)')
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
