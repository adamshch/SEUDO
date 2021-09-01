function [theEval,transArea] = evaluateOneTimeCourse(tcOrig,tcNew,varargin)


p = inputParser;
p.addParamValue('method','area');
p.addParamValue('threshSigTrans',0.1); % threshold for being considered a significant part of the transient (dF/F units)
p.addParamValue('threshTall',0.2); % how far (fraction) the new time course can fall and still be considered tall enough
p.addParamValue('minDuration',0.5); % fraction of original duration time course can be shortened to still count as long enough
parse(p,varargin{:});

params = p.Results;


transArea = sum(tcOrig(tcOrig>0));

switch params.method
    case 'area' % fraction kept
        
        tcOrig(tcOrig<0) = 0;
        tcNew(tcNew<0) = 0;
        
        theEval = sum(tcNew)/sum(tcOrig);
        
    case 'tall span' % fraction of points above X% of original
        
        thresh = params.threshSigTrans;
        fallFraction = params.threshTall;
        
        ptsToCheck = tcOrig > thresh;
        
        if any(ptsToCheck)
            theEval = sum(tcNew(ptsToCheck)>tcOrig(ptsToCheck)*fallFraction)/sum(ptsToCheck);
        else
            theEval = nan;
        end
        
        
    case 'tall span threshold' % duration greater than 
        
        thresh = params.threshSigTrans;
        minDuration = params.minDuration;
        fallFraction = params.threshTall;
        
        ptsToCheck = tcOrig > thresh;
        
        if any(ptsToCheck)
            theEval = sum(tcNew(ptsToCheck)>tcOrig(ptsToCheck)*fallFraction) >= sum(ptsToCheck)*minDuration;
        else
            theEval = nan;
        end
        
    otherwise
        error('method %s not recognized',params.method)
end
