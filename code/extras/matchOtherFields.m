function [s1,s2] = matchOtherFields(s1, s2)                                                                                                    

% s1 = matchOtherFields(s1, s2)
%
%
% 2019 - Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Start with default parameters

if numel(s2) > 1
    %s2_new = [];
    for ll = 1:numel(s2)
        [s1,s2_new(ll)] = matchOtherFields(s1, s2(ll));
    end
    s2 = s2_new;
else
    f2 = fieldnames(s2);
    for ll = 1:numel(f2)
        if ~isfield(s1,f2{ll})
            s1.(f2{ll}) = [];
        end
    end
    f1 = fieldnames(s1);
    for ll = 1:numel(f1)
        if ~isfield(s2,f1{ll})
            [s2.(f1{ll})] = deal(cell(1,numel(s2)));
        end
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
