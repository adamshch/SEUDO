function varargout = plotCorrDetectROC(se,varargin)

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get bounds to sweep over

min_corr = 0;
max_corr = 1;
for ll = 1:numel(se)
    for cc = 1:se{ll}.nCells
        min_corr = min(min_corr, ...
                 min(se{ll}.tcDefault.transientInfo(cc).corrWithProfile));
        max_corr = max(max_corr, ...
                 max(se{ll}.tcDefault.transientInfo(cc).corrWithProfile));
    end
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
    for ll = 1:numel(se)
        CA = se{ll}.autoClassifyTransients('default', 'saveResults', ...
                     false,'plotComparison',false, 'method', 'corr', ...
                      'ignoreZeros',false, 'corrThresh', thresh_vals(kk)); % Calculate the threshold classification
        for cc = 1:se{ll}.nCells
            classHuman = cat(1,classHuman,...
                       se{ll}.tcDefault.transientInfo(cc).classification); % Store the human classifications
            classAuto  = cat(1,classAuto,CA(cc).classification);           % Get the new automatic classifications
        end
    end
    p_tp(kk) = sum((classHuman== seudo.valTrue)  & (classAuto==seudo.valTrue)) / sum(classHuman== seudo.valTrue);  	% Count the true positives
    p_fp(kk) = sum((classHuman== seudo.valFalse) & (classAuto==seudo.valTrue)) / sum(classHuman== seudo.valFalse); 	% Count the false positives
    p_am(kk) = sum((classHuman== seudo.valMix)   & (classAuto==seudo.valTrue)) / sum(classHuman== seudo.valMix);    % Count the mixed positives
    p_un(kk) = sum( isnan(classHuman)            & (classAuto==seudo.valTrue)) / sum(isnan(classHuman));            % Count the unclassified positives
    fprintf('.')
end
fprintf('\n')

h = figure;
subplot(2,3,[1,5]), plot(p_fp,p_tp,'g')
xlabel('P(false positive)')
ylabel('P(true positive)')
box off
axis square
set(gca,'TickDir','out','XTick',0:0.2:1,'YTick',0:0.2:1)


subplot(2,3,3), plot(thresh_vals,p_tp,'b')
xlabel('Correlation threshold')
yyaxis left
ylabel('P(true positive)')
yyaxis right
subplot(2,3,3), plot(thresh_vals,p_fp,'r')
xlabel('Correlation threshold')
ylabel('P(false positive)')
box off
set(gca,'TickDir','out','YTick',0:0.5:1)

subplot(2,3,6), plot(thresh_vals,p_am,'g')
xlabel('Correlation threshold')
yyaxis left
ylabel('P(mixed positive)')
yyaxis right
subplot(2,3,6), plot(thresh_vals,p_un,'color', [0.5,0.5,0.5])
xlabel('Correlation threshold')
ylabel('P(unclassified positive)')
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