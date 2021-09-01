
cT = 0.5;

[mc, mr, countPerCell] = calcFullAutoClass(se_cell,'method', 'corr','ignoreZeros',false, 'weightedTrans', true,'resRatioMean',false, 'corrThresh', cT);
plotFaLCon(mc, mr, 'counts', sum(countPerCell,2), 'numTrue', countPerCell(:,1), 'numClass', sum(countPerCell(:,1:3),2));


frac_tm = [];
for ll = 1:numel(se_cell)
    for cc = 1:se_cell{ll}.nCells
        se_class = se_cell{ll}.tcDefault.transientInfo(cc).classification;
        se_corr  = se_cell{ll}.tcDefault.transientInfo(cc).corrWithProfile;
        frac_tm  = cat(1, frac_tm, sum((se_corr>=cT)&((se_class==1)|(se_class==0)))/numel(se_class));
%         frac_rz  = cat(1, frac_rz, se_cell{ll}.tcDefault.transientInfo(cc). );
    end
end

figure
scatter(frac_tm, mr)
title(sprintf('Corr = %f', cT))
xlabel('Fraction true or mixed')
ylabel('Contam Severity')
% clear se_class se_corr


%%

% Possible_UL = [1078];
% Possible_UR = [228];
% Possible_LL = [865];
% Possible_LR = [454];

% ex_cells = [1078, 103, 1004, 916];                                          % Example cells in: UL,UR,LL, LR,

ex_cells = [];
% ex_cells = [973, 1237, 1481, 880];                                          % Example cells in: UL,UR,LL, LR,
% ex_cells = [];
% plotFaLConRandExamples(se_cell, mc, mr, 'ex_cells' , [], 'plotClassOnly', true, 'figNum', 2)
plotFaLConRandExamples(se_cell, mc, mr, 'ex_cells' , ex_cells, 'plotClassOnly', false, 'figNum', 2)

