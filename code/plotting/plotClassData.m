function varargout = plotClassData(se, varargin)

% varargout = plotClassData(se, varargin)
%
% Function to plot histograms showing the amplitudes and correlations of
% differently classified transients (via SEUDO objects)
%
% 2018 - Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input parsing

if ~iscell(se)
    TMP   = se;
    se    = cell(1);
    se{1} = TMP;
    clear TMP
end

if nargin > 1
    n_bins = varargin{1};
else
    n_bins = 20;
end

if nargin > 2
    cellsToPlot = varargin{2};
else
    cellsToPlot = cell(size(se));
end

if ~iscell(cellsToPlot)
    TMP            = cellsToPlot;
    cellsToPlot    = cell(1);
    cellsToPlot{1} = TMP;
    clear TMP
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get true/false amplitudes

amp_true  = [];
amp_false = [];
amp_mixed = [];
amp_uncls = [];

res_true  = [];
res_false = [];
res_mixed = [];
res_uncls = [];


cor_true  = [];
cor_false = [];
cor_mixed = [];
cor_uncls = [];


for ll = 1:numel(se)
    if isempty(cellsToPlot{ll})
        cellsToPlot{ll} = 1:size(se{ll}.tcDefault.tc,2);
    end
%     for cc = 1:size(se{ll}.tcDefault.tc,2)
    for cc = cellsToPlot{ll}
        
        evts     = se{ll}.tcDefault.transientInfo(cc).times;
        classVal = se{ll}.tcDefault.transientInfo(cc).classification;
        tc_corrs = se{ll}.tcDefault.transientInfo(cc).corrWithProfile;
        r        = calcResRatios(se{ll}.profiles, ...
                                  se{ll}.tcDefault.transientInfo(cc), cc);
        tc_amps  = zeros(size(evts,1),1);                                  % Initialize array of amplitudes
        for kk = 1:size(evts,1)
            tc_amps(kk) = max(se{ll}.tcDefault.tc(evts(kk,1):evts(kk,2),cc));
        end
        
        
        amp_true  = cat(1,amp_true, vec(tc_amps(classVal==seudo.valTrue)));  	% Separate out amplitudes for true transients...
        amp_false = cat(1,amp_false,vec(tc_amps(classVal==seudo.valFalse)));   	%    ... false transients ...
        amp_mixed = cat(1,amp_mixed,vec(tc_amps(classVal==seudo.valUnc )));    	%    ... mixed transients ...
        amp_uncls = cat(1,amp_uncls,vec(tc_amps(isnan(classVal))));             %    ... unclassified transients.        
        
        res_true  = cat(1,res_true, vec(r.resRatios(classVal==seudo.valTrue )));% Separate out amplitudes for true transients...
        res_false = cat(1,res_false,vec(r.resRatios(classVal==seudo.valFalse)));%    ... false transients ...
        res_mixed = cat(1,res_mixed,vec(r.resRatios(classVal==seudo.valUnc ))); %    ... mixed transients ...
        res_uncls = cat(1,res_uncls,vec(r.resRatios(isnan(classVal))));         %    ... unclassified transients.
        
        cor_true  = cat(1,cor_true, vec(tc_corrs(classVal==seudo.valTrue ))); 	% Separate out correlation values for true transients...
        cor_false = cat(1,cor_false,vec(tc_corrs(classVal==seudo.valFalse)));  	%    ... false transients ...
        cor_mixed = cat(1,cor_mixed,vec(tc_corrs(classVal==seudo.valUnc)));    	%    ... mixed transients ...
        cor_uncls = cat(1,cor_uncls,vec(tc_corrs(isnan(classVal))));            %    ... unclassified transients.
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Perform plotting

h(1) = figure;
subplot(4,3,1), histogram(amp_true, n_bins) % histfit(amp_true, n_bins,'kernel')
box off
set(gca, 'TickDir','out')
ylabel('Frequency')
xlabel('True amplitudes')
subplot(4,3,4), histogram(amp_false, n_bins)
box off
set(gca, 'TickDir','out')
ylabel('Frequency')
xlabel('False amplitudes')
subplot(4,3,7), histogram(amp_mixed, n_bins)
box off
set(gca, 'TickDir','out')
ylabel('Frequency')
xlabel('mixed amplitudes')
subplot(4,3,10), histogram(amp_uncls, n_bins)
box off
set(gca, 'TickDir','out')
ylabel('Frequency')
xlabel('Unclassified amplitudes')

subplot(4,3,2), histogram(cor_true, n_bins)
box off
set(gca, 'TickDir','out')
ylabel('Frequency')
xlabel('True correlations')
subplot(4,3,5), histogram(cor_false, n_bins)
box off
set(gca, 'TickDir','out')
ylabel('Frequency')
xlabel('False correlations')
subplot(4,3,8), histogram(cor_mixed, n_bins)
box off
set(gca, 'TickDir','out')
ylabel('Frequency')
xlabel('mixed correlations')
subplot(4,3,11), histogram(cor_uncls, n_bins)
box off
set(gca, 'TickDir','out')
ylabel('Frequency')
xlabel('Unclassified correlations')


subplot(4,3,3), histogram(res_true, n_bins)
box off
set(gca, 'TickDir','out')
ylabel('Frequency')
xlabel('True residual ratio')
subplot(4,3,6), histogram(res_false, n_bins)
box off
set(gca, 'TickDir','out')
ylabel('Frequency')
xlabel('False residual ratio')
subplot(4,3,9), histogram(res_mixed, n_bins)
box off
set(gca, 'TickDir','out')
ylabel('Frequency')
xlabel('mixed residual ratio')
subplot(4,3,12), histogram(res_uncls, n_bins)
box off
set(gca, 'TickDir','out')
ylabel('Frequency')
xlabel('Unclassified residual ratio')

set(gcf,'color',[1,1,1])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cumulative histogram plots

HA1 = fitdist(amp_true,  'Kernel', 'Kernel', 'epanechnikov', 'Support', 'positive');
HA2 = fitdist(amp_false, 'Kernel', 'Kernel', 'epanechnikov', 'Support', 'positive');
HA3 = fitdist(amp_mixed, 'Kernel', 'Kernel', 'epanechnikov', 'Support', 'positive');
HA4 = fitdist(amp_uncls, 'Kernel', 'Kernel', 'epanechnikov', 'Support', 'positive');

HC1 = fitdist(cor_true,  'Kernel', 'Kernel', 'epanechnikov');
HC2 = fitdist(cor_false, 'Kernel', 'Kernel', 'epanechnikov');
HC3 = fitdist(cor_mixed, 'Kernel', 'Kernel', 'epanechnikov');
HC4 = fitdist(cor_uncls, 'Kernel', 'Kernel', 'epanechnikov');

HR1 = fitdist(res_true,  'Kernel', 'Kernel', 'epanechnikov', 'Support', 'positive');
HR2 = fitdist(res_false, 'Kernel', 'Kernel', 'epanechnikov', 'Support', 'positive');
HR3 = fitdist(res_mixed, 'Kernel', 'Kernel', 'epanechnikov', 'Support', 'positive');
HR4 = fitdist(res_uncls, 'Kernel', 'Kernel', 'epanechnikov', 'Support', 'positive');

h(2) = figure;
x_vals = linspace(0*min([amp_true(:); amp_false(:); amp_mixed(:); amp_uncls(:)]), ...
                  1.1*max([amp_true(:); amp_false(:); amp_mixed(:); amp_uncls(:)]), 100);
subplot(2,3,1), plot(x_vals, pdf(HA1,x_vals), 'b') 
subplot(2,3,1), hold on
subplot(2,3,1), plot(x_vals, pdf(HA2,x_vals), 'r') 
subplot(2,3,1), plot(x_vals, pdf(HA3,x_vals), 'g') 
subplot(2,3,1), plot(x_vals, pdf(HA4,x_vals), 'color', 0.3*[1,1,1]) 
subplot(2,3,1), hold off
box off
set(gca, 'TickDir','out')
legend('True', 'False', 'Mixed','Artifact')
ylabel('Probability')
xlabel('Amplitude')

subplot(2,3,4), plot(x_vals, cdf(HA1,x_vals), 'b') 
subplot(2,3,4), hold on
subplot(2,3,4), plot(x_vals, cdf(HA2,x_vals), 'r') 
subplot(2,3,4), plot(x_vals, cdf(HA3,x_vals), 'g') 
subplot(2,3,4), plot(x_vals, cdf(HA4,x_vals), 'color', 0.3*[1,1,1]) 
subplot(2,3,4), hold off
box off
set(gca, 'TickDir','out')
legend('True', 'False', 'Mixed','Artifact')
ylabel('Cumulative probability')
xlabel('Amplitude')


x_vals = linspace(1.1*min([cor_true(:); cor_false(:); cor_mixed(:); cor_uncls(:)]), ...
                  max([cor_true(:); cor_false(:); cor_mixed(:); cor_uncls(:)]), 100);
subplot(2,3,2), plot(x_vals, pdf(HC1,x_vals), 'b') 
subplot(2,3,2), hold on
subplot(2,3,2), plot(x_vals, pdf(HC2,x_vals), 'r') 
subplot(2,3,2), plot(x_vals, pdf(HC3,x_vals), 'g') 
subplot(2,3,2), plot(x_vals, pdf(HC4,x_vals), 'color', 0.3*[1,1,1]) 
subplot(2,3,2), hold off
box off
set(gca, 'TickDir','out')
legend('True', 'False', 'Mixed','Artifact')
ylabel('Probability')
xlabel('Correlation')

subplot(2,3,5), plot(x_vals, cdf(HC1,x_vals), 'b') 
subplot(2,3,5), hold on
subplot(2,3,5), plot(x_vals, cdf(HC2,x_vals), 'r') 
subplot(2,3,5), plot(x_vals, cdf(HC3,x_vals), 'g') 
subplot(2,3,5), plot(x_vals, cdf(HC4,x_vals), 'color', 0.3*[1,1,1]) 
subplot(2,3,5), hold off
box off
set(gca, 'TickDir','out')
legend('True', 'False', 'Mixed','Artifact')
ylabel('Cumulative probability')
xlabel('Correlation')


x_vals = linspace(0*min([res_true(:); res_false(:); res_mixed(:); res_uncls(:)]), ...
                  1.1*max([res_true(:); res_false(:); res_mixed(:); res_uncls(:)]), 100);
subplot(2,3,3), plot(x_vals, pdf(HR1,x_vals), 'b') 
subplot(2,3,3), hold on
subplot(2,3,3), plot(x_vals, pdf(HR2,x_vals), 'r') 
subplot(2,3,3), plot(x_vals, pdf(HR3,x_vals), 'g') 
subplot(2,3,3), plot(x_vals, pdf(HR4,x_vals), 'color', 0.3*[1,1,1]) 
subplot(2,3,3), hold off
box off
set(gca, 'TickDir','out')
legend('True', 'False', 'Mixed','Artifact')
ylabel('Probability')
xlabel('Residual ratio')

subplot(2,3,6), plot(x_vals, cdf(HR1,x_vals), 'b') 
subplot(2,3,6), hold on
subplot(2,3,6), plot(x_vals, cdf(HR2,x_vals), 'r') 
subplot(2,3,6), plot(x_vals, cdf(HR3,x_vals), 'g') 
subplot(2,3,6), plot(x_vals, cdf(HR4,x_vals), 'color', 0.3*[1,1,1]) 
subplot(2,3,6), hold off
box off
set(gca, 'TickDir','out')
legend('True', 'False', 'Mixed','Artifact')
ylabel('Cumulative probability')
xlabel('Residual ratio')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output parsing

if nargout == 1
    varargout{1} = h;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
