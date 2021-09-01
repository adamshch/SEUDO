

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get pointer to file with all the time-traces

% Original file at: /jukebox/tank/jlgauthi/SEUDO/J122_2015-11-16_L01_c01_tc.mat
%                   /jukebox/tank/jlgauthi/SEUDO/J122_2015-11-16_L01_c01_tc.mat
m = matfile('/mnt/bucket/adam/SEUDOdata/J122_2015-11-16_L01_c01_tc.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Extract the time traces

tcLSQ = single(m.tcAllLSQ);                                                % Least squares time-traces
tcCor = single(m.tcAllCorr);                                               % Correlation corrected time-traces
tcSEU = single(m.tcAllSEUDO);                                              % SEUDO-corrected time-traces

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate covariance, correlation matrices and their spectra

cov.tcCor = tcCor.'*tcCor - mean(tcCor,1).'*mean(tcCor,1);
cov.tcLSQ = tcLSQ.'*tcLSQ - mean(tcLSQ,1).'*mean(tcLSQ,1);
cov.tcSEU = tcSEU.'*tcSEU - mean(tcSEU,1).'*mean(tcSEU,1);

corr.tcCor = diag(sqrt(1./diag(cov.tcCor)))*cov.tcCor*diag(sqrt(1./diag(cov.tcCor)));
corr.tcLSQ = diag(sqrt(1./diag(cov.tcLSQ)))*cov.tcLSQ*diag(sqrt(1./diag(cov.tcLSQ)));
corr.tcSEU = diag(sqrt(1./diag(cov.tcSEU)))*cov.tcSEU*diag(sqrt(1./diag(cov.tcSEU)));

eigMat = [abs(svd(corr.tcCor)),abs(svd(corr.tcLSQ)),abs(svd(corr.tcSEU))];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Clustering

clust.num = 5;
clust.LSQ = kmeans(corr.tcLSQ.',clust.num,'MaxIter',5000,'Replicates',500); 
clust.Cor = kmeans(corr.tcCor.',clust.num,'MaxIter',5000,'Replicates',500); 
clust.SEU = kmeans(corr.tcSEU.',clust.num,'MaxIter',5000,'Replicates',500); 

%%
clust.LSQ = sortClusters(clust.LSQ);
[clust.LSQ,clust.Cor,clust.confLC] = relableMatches(clust.LSQ,clust.Cor);
[clust.LSQ,clust.SEU,clust.confLS] = relableMatches(clust.LSQ,clust.SEU);
clust.confCS                       = confusionmat(clust.Cor,clust.SEU);

%%

figure(1);
subplot(3,3,1), cla; histogram(clust.LSQ)
box off; title('Least-squares'); xlabel('Label'); ylabel('Frequency')
subplot(3,3,2), cla; histogram(clust.Cor)
box off; title('Correlation'); xlabel('Label'); ylabel('Frequency')
subplot(3,3,3), cla; histogram(clust.SEU)
box off; title('SEUDO'); xlabel('Label'); ylabel('Frequency')

[~,ix.LSQ] = sort(clust.LSQ,'ascend');
[~,ix.COR] = sort(clust.Cor,'ascend');
[~,ix.SEU] = sort(clust.SEU,'ascend');

subplot(3,3,4), cla; imagesc(corr.tcLSQ(ix.LSQ,ix.LSQ).*(1-eye(numel(ix.LSQ))))
axis image; box off; title('Least-squares'); xlabel('Neurons'); ylabel('Neurons')
subplot(3,3,5), cla; imagesc(corr.tcCor(ix.COR,ix.COR).*(1-eye(numel(ix.LSQ))))
axis image; box off; title('Correlation'); xlabel('Neurons'); ylabel('Neurons')
subplot(3,3,6), cla; imagesc(corr.tcSEU(ix.SEU,ix.SEU).*(1-eye(numel(ix.LSQ))))
axis image; box off; title('SEUDO'); xlabel('Neurons'); ylabel('Neurons')

subplot(3,3,7), confusionchart(clust.confLC,'Normalization','row-normalized')
title('LS/Corr confusion'); xlabel('Least-squares'); ylabel('Correlation')
subplot(3,3,8), confusionchart(clust.confLS,'Normalization','row-normalized')
title('LS/SEUDO confusion'); xlabel('Least-squares'); ylabel('SEUDO')
subplot(3,3,9), confusionchart(clust.confCS,'Normalization','row-normalized')
title('Corr/SEUDO confusion'); xlabel('Correlation'); ylabel('SEUDO')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot the correlation matrices nad the singular value spectra

figure(2); plot(bsxfun(@times,cumsum(eigMat),1./sum(eigMat,1)))
set(gca,'Yscale','lin','Xscale','log'); box off
title('Singular value spectrum')
legend('LSQ','Correlation','SEUDO')
xlabel('Singular value number')
ylabel('Log-cumulative variance explained')

figure(3); plot(eigMat)
set(gca,'Yscale','lin','Xscale','log'); box off
title('Singular value spectrum')
legend('LSQ','Correlation','SEUDO')
xlabel('Singular value number')
ylabel('Log-cumulative variance explained')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Artifact trace plotting

isArt     = m.isArtifact;
idxArt    = find(isArt==1);
tiLSQ     = seudo.identifyTransients(tcLSQ);
tiSEU     = seudo.identifyTransients(tcSEU);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 

idxToPlot = 11;
ratioArt  = sum(tiSEU(:,isArt))./sum(tiLSQ(:,isArt));
ratioNot  = sum(tiSEU(:,~isArt))./sum(tiLSQ(:,~isArt));
figure(4); 
subplot(2,1,1), cla; subplot(2,1,1), hold on
histogram(ratioArt,10)
plot(nanmean(ratioArt)*[1,1],[0,40], '--r')
set(gca,'Xlim',[0,1.1]); box off
title('Artifact cells')
ylabel('Number of profiles')
xlabel('(active SEUDO frames)/(active LSQ frames)')
subplot(2,1,2), cla; subplot(2,1,2), hold on
histogram(ratioNot,10)
plot(nanmean(ratioNot)*[1,1],[0,300], '--r')
set(gca,'Xlim',[0,1.1]); box off
title('Non-artifact cells')
ylabel('Number of profiles')
xlabel('(active SEUDO frames)/(active LSQ frames)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%