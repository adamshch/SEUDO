function [m_cfrac, m_rfrac, countPerCell, varargout] = calcFullAutoClass(se, varargin)

% [m_cfrac, m_rfrac, countPerCell] = calcFullAutoClass(se_cell, varargin)
%
% 
%
% 2018 - Adam Charles 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 

if ~iscell(se)                                                             % This block makes sure that se is a cell
    TMP   = se;
    se    = cell(1);
    se{1} = TMP;
    clear TMP
end

m_cfrac       = [];
m_rfrac       = [];
m_rfrac2      = [];
m_rfrac3      = [];

countPerCell  = [];
for kk = 1:numel(se)
    fprintf('Loading block %d...',kk)
    CA = se{kk}.autoClassifyTransients('default','saveResults',...
                              false, 'plotComparison',false,varargin{:});  % Claculate automatic classification
    m_cfrac = cat(1,m_cfrac,[CA.cfrac].');
    m_rfrac = cat(1,m_rfrac,[CA.rfrac].');
    m_rfrac2 = cat(1,m_rfrac2,[CA.rfrac2].');
    m_rfrac3 = cat(1,m_rfrac3,[CA.rfrac3].');
    
    for cc = 1:se{kk}.nCells
        C             = se{kk}.tcDefault.transientInfo(cc).classification; % Pull out classification data
        countPerCell  = cat(1, countPerCell, ...
                        [sum(C==1) sum(C==0) sum(C==-1), sum(isnan(C))]);  % Collect different classification counts
    end
    clear CA
    fprintf('done.\n')
end

if nargout > 3
    varargout{1} = m_rfrac2;
end
if nargout > 4
    varargout{2} = m_rfrac3;
end

end

%     CA        = se_cell{kk}.autoClassifyTransients('default','saveResults',false,...
%                       'plotComparison',false,'radius',[3 7],'alpha',0.01,'tauBounds',[0.05,0.1]); % automatically classify transients
%     CA = se_cell{kk}.autoClassifyTransients('default','saveResults',...
%                      false, 'plotComparison',false,'method', 'corr',...
%                                 'ignoreZeros',false, 'corrThresh', 0.5);  % Claculate automatic classification


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%