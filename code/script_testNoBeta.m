%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

% m1 = matfile('/home/adam/BigStorage/Data/SEUDO/SEUDO_J122_2015-11-16_L01_c01_b27.mat');
% m2 = matfile('/home/adam/BigStorage/Data/SEUDO/J122_2015-11-16_L01_c01_b27_class_ASC2.mat');
% m3 = matfile('/home/adam/BigStorage/Data/SEUDO/caiman_SEUDO_J122_2015-11-16_L01_c01_b27.mat');
seudoData = load('/home/adam/BigStorage/Data/SEUDO/SEUDO_J122_2015-11-16_L01_c01_b27_ForClass.mat');
% seudo(seudoData.F_movie,seudoData.cellFile.rois,'timeCourses',...
%                seudoData.cellFile.timeCourses','tcDefault',fileclass{kk});
% fileclass = '/home/adam/BigStorage/Data/SEUDO/J122_2015-11-16_L01_c01_b27_class_ASC2.mat'; % J122_2015-11-16_L01_c01_b27_ACclass.mat
fileclass = '/home/adam/BigStorage/Data/SEUDO/J122_2015-11-16_L01_c01_b27_ACclass.mat';
% caiman_out = m3.caiman_out;
% TMP = m1.cellFile;

seudoData = load('/home/adam/BigStorage/Data/SEUDO/Suite2pClass/SEUDO_Suite2pData.mat');

lambdaList = [10, 100, 300, 500, 600, 700, 800, 900, 1000, 2000, 3000, 1e4, 1e5];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

% se = seudo(M,P,'timeCourses',T,'name','demo dataset','tcDefault',tDefaultDemo);
% se = seudo(seudoData.F_movie,seudoData.cellFile.rois,'timeCourses',...
%                seudoData.cellFile.timeCourses','name','demo dataset','tcDefault',fileclass);
se = seudo(double(seudoData.M),double(seudoData.A),'timeCourses',...
               double(seudoData.T)','name','demo dataset','tcDefault',seudoData.se.tcDefault);
           
% se = seudoData.se;
% se.tcSeudo = [];
           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

for ll = lambdaList
    seudoParams = {'dsTime',3,'sigma2',0.0020,'lambdaBlob',ll,'blobRadius',3,'padSpace',5,'saveBlobTimeCourse',1};
    se.parallelSEUDO([], seudoParams{:})
end

results1 = se.computeSeudoROC;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

% se2 = seudo(M,P,'timeCourses',T,'name','demo dataset','tcDefault',tDefaultDemo);
% se2 = seudo(m1.F_stab,TMP.rois,'timeCourses',TMP.timeCourses.','name','demo dataset','tcDefault',saveName);
se2 = seudo(double(seudoData.M),double(seudoData.A),'timeCourses',...
               double(seudoData.T)','name','demo dataset','tcDefault',seudoData.se.tcDefault);

% se2 = seudoData.se;
se2.tcSeudo = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

for ll = lambdaList
    seudoParams = {'dsTime',3,'sigma2',0.0020,'lambdaBlob',ll,'blobRadius',3,'padSpace',5,'saveBlobTimeCourse',1};
    se2.parallelSEUDO([], seudoParams{:})
end

results2 = se2.computeSeudoROC;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

results1 = se.computeSeudoROC;
results2 = se2.computeSeudoROC;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

subplot(1,2,1), plot([results1.meanAreaFalse], [results1.meanAreaTrue], '-b', [results2.meanAreaFalse], [results2.meanAreaTrue], '--r')
box off
xlabel('False transient area kept')
ylabel('True transient area kept')
legend('With beta', 'No beta')
set(gca, 'XLim', [0,1], 'YLim', [0, 1])

tt1 = zeros(numel(results1),1);
ft1 = zeros(numel(results1),1);
tt2 = zeros(numel(results1),1);
ft2 = zeros(numel(results1),1);
for ll = 1:numel(results1)
    ft1(ll) = mean(results1(ll).falseTransDetected);
    tt1(ll) = mean(results1(ll).trueTransDetected);
    ft2(ll) = mean(results2(ll).falseTransDetected);
    tt2(ll) = mean(results2(ll).trueTransDetected);
end

subplot(1,2,2), plot(ft1(1:end), tt1(1:end), '-b', ft2(1:end), tt2(1:end), '--r')
set(gca, 'XLim', [0,1], 'YLim', [0, 1])
xlabel('False transients kept')
ylabel('True transients kept')
legend('With beta', 'No beta')
box off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

ft1b = [1; ft1(:); 0];
ft2b = [1; ft2(:); 0];
tt1b = [1; tt1(:); 0];
tt2b = [1; tt2(:); 0];

AUC1wb = sum(((ft1b(1:end-1) - ft1b(2:end)).*0.5.*(tt1b(2:end)+tt1b(1:end-1)) ) );
AUC1nb = sum(((ft2b(1:end-1) - ft2b(2:end)).*0.5.*(tt2b(2:end)+tt2b(1:end-1)) ) );

ft1b = [1; vec([results1.meanAreaFalse]); 0];
ft2b = [1; vec([results2.meanAreaFalse]); 0];
tt1b = [1; vec([results1.meanAreaTrue]); 0];
tt2b = [1; vec([results2.meanAreaTrue]); 0];

AUC2wb = sum(((ft1b(1:end-1) - ft1b(2:end)).*0.5.*(tt1b(2:end)+tt1b(1:end-1)) ) );
AUC2nb = sum(((ft2b(1:end-1) - ft2b(2:end)).*0.5.*(tt2b(2:end)+tt2b(1:end-1)) ) );

fprintf('AUC values:\n')
fprintf('beta>0:  (pct) = %f, (area) = %f\n', AUC1wb, AUC2wb)
fprintf('beta=0:  (pct) = %f, (area) = %f\n', AUC1nb, AUC2nb)
fprintf('Delta:  (pct) = %f, (area) = %f\n', AUC1wb - AUC1nb, AUC2wb - AUC2nb)

clear ft1b ft2b tt1b tt2b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

