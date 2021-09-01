function results = computeSeudoROC(se,varargin)
% seudoObject.computeSeudoROC(...)
%
% evaluate SEUDO performance on transients classified as true or false
%
%
% EXAMPLES
%
%   se.computeSeudoROC
%   se.computeSeudoROC('whichSeudo',5:9,'whichDefault',1)
%   pts = se.computeSeudoROC
%   


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% parse inputs
p = inputParser;
p.addParameter('whichCells',[]);
p.addParameter('whichSeudo',1:length(se.tcSeudo));
p.addParameter('whichDefault',1);
p.addParameter('plot',true);
parse(p,varargin{:});


params = p.Results;
clear p


% parse user input to select time course structs
structSpecOrig = pickTcStruct(se,{'default',params.whichDefault});
%structSpecSeudo = pickTcStruct(se,{'seudo',params.whichSeudo});


% ensure transients are classified in original time course struct
if ~isfield(se.(structSpecOrig{1})(structSpecOrig{2}),'transientInfo') || ...
        isempty(se.(structSpecOrig{1})(structSpecOrig{2}).transientInfo)
    error('must compute transient info and classify transients first')
end

% ensure transient info is computed for seudo time courses
% if ~isfield(se.(structSpecSeudo{1})(structSpecSeudo{2}),'transientInfo') || ...
%         isempty(se.(structSpecSeudo{1})(structSpecSeudo{2}).transientInfo)
%     se.computeTransientInfo({'seudo',params.whichSeudo})
% end


% loop over SEUDO parameter sets
for ss = 1:length(params.whichSeudo)
    
    % select seudo time course struct
    structSpecSeudo = pickTcStruct(se,{'seudo',params.whichSeudo(ss)});
    
    
    % for every transient, get classification, area, and whether detected
    allTransClass  = [];
    allTransArea   = [];
    allTransDetect = [];
    
    
    % loop over cells...
    for cc = 1:se.nCells
        
        % get transient info for this cell from original and seudo
        tiOrig = se.(structSpecOrig{1})(structSpecOrig{2}).transientInfo(cc);
        
        % get time courses
        tcOrig  = se.(structSpecOrig{1})(structSpecOrig{2}).tc(:,cc);
        tcLSQ   = se.(structSpecSeudo{1})(structSpecSeudo{2}).tcLSQ(:,cc);
        tcSeudo = se.(structSpecSeudo{1})(structSpecSeudo{2}).tc(:,cc);
        
        % make time course which is SEUDO during transients but otherwise same as the original
        %   (this is needed because we want to run the same transient classification on the
        %   seudo data as the original, and this is based on the noise level from non-transient frames)
        
        
        % find scale factor between computing seudo time course and original
        %   (needed in case the absolute values of the movie changed between
        %   computing original time courses and computing seudo)
        analyzedFrames = ~isnan(tcLSQ);
%         scaleFactor    = tcLSQ(analyzedFrames) \ tcOrig(analyzedFrames);
        scaleFactor    = median(tcOrig(analyzedFrames)./tcLSQ(analyzedFrames));
        
        % fill in scaled-up version of SEUDO time course
        tcFranken = tcOrig;
        tcFranken(analyzedFrames) = 1.3*tcSeudo(analyzedFrames)*scaleFactor;
        
        % identify transients in SEUDO time course using same parameters as original
        seudoTrans = se.identifyTransients(tcFranken,...
                          'threshFcn',   tiOrig.params.threshFcn,...
                          'blurRadius',  tiOrig.params.blurRadius,...
                          'minDuration', tiOrig.params.minDuration);
        
        
        % plot to debug
        if 0
            figure(101);clf;
            subplot(1,2,1);plot(tcOrig);hold on;plot(tcLSQ,'mo-')
            subplot(1,2,2);plot(tcOrig,tcLSQ,'.');%hold on;plot(xlim,xlim,'k-')
            pause
        end
        
        
        % identify area for each transient from original time course
        for tt = 1:length(tiOrig.classification)
            
            % note frames of this transient
            whichFrames = tiOrig.times(tt,1):tiOrig.times(tt,2);
            
            % get original and seudo time course
            transOrig  = tcLSQ(whichFrames);
            transSeudo = tcSeudo(whichFrames);
            
            % rectify
            transOrig(transOrig<0)   = 0;
            transSeudo(transSeudo<0) = 0;
            
            % compute active fraction
            tcArea  = sum(transSeudo) / sum(transOrig);
            
            % get classification
            tcClass = tiOrig.classification(tt);
            
            % note whether transient was detected
            tcDetected = any(seudoTrans(whichFrames));
            
            % store results
            allTransClass  = [allTransClass;  tcClass];
            allTransArea   = [allTransArea;   tcArea];
            allTransDetect = [allTransDetect; tcDetected];
            
        end
    end
    
    
    % identify true and false transients (from original time course)
    trueTrans  = allTransClass == seudo.getClassValue('true');
    falseTrans = allTransClass == seudo.getClassValue('false');
    
    
    % get mean areas
    meanAreaTrue  = nanmean(allTransArea(trueTrans));
    meanAreaFalse = nanmean(allTransArea(falseTrans));
    
    
    % package output
    theseResults = struct;
    theseResults.trueTransAreas     = allTransArea(trueTrans);
    theseResults.falseTransAreas    = allTransArea(falseTrans);
    theseResults.meanAreaTrue       = meanAreaTrue;
    theseResults.meanAreaFalse      = meanAreaFalse;
    theseResults.trueTransDetected  = allTransDetect(trueTrans);
    theseResults.falseTransDetected = allTransDetect(falseTrans);
    
    % append to list of results
    results(ss) = theseResults;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot

if params.plot
    
    figure%(1760);clf
    set(gcf,'name',sprintf('%s ROC plot',se.name))
    
    % area preserved
    subplot(1,2,1)
    x = [results.meanAreaFalse];
    y = [results.meanAreaTrue];
    plot(x,y,'.k','markersize',15)
    for ss = 1:length(x)
        text(x(ss),y(ss),['  ' num2str(ss)],'horizontalalignment','left')
    end
    title('area of transients preserved')
    
    % fraction of transients detected
    subplot(1,2,2)
    ftd = cellfun(@mean,{results.falseTransDetected});
    ttd = cellfun(@mean,{results.trueTransDetected});
    %plot(mean(theseResults.falseTransDetected),mean(theseResults.trueTransDetected),'k.','markersize',15)
    plot(ftd,ttd,'k.','markersize',15)
    for ss = 1:length(ftd)
        text(ftd(ss),ttd(ss),['  ' num2str(ss)],'horizontalalignment','left')
    end
    title('fraction of transients detected')
    
    % make pretty
    for pp = 1:2
        subplot(1,2,pp)
        
        set(gca,'xtick',0:.2:1,'ytick',0:.2:1)
        grid on
        axis image
        xlim([0 1])
        ylim([0 1])
        xlabel('false');ylabel('true')
    end
    
end

% don't return if not requested
if nargout < 1, clear results, end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%