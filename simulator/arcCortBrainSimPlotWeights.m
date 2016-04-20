% This script takes a fiber group and signal and outputs plots of the
% original fiber weights and LiFE weights

% Sam Faber 2016/03/21

%% Set up paths and build the file names for the diffusion data, bvecs, bvals, and fiber group. 
subj      = '105115';
rootDir   = fullfile('/N/dc2/projects/lifebid/HCP/Sam/',subj);
saveDir   = fullfile('/N/dc2/projects/lifebid/HCP/Sam/matlab_code/wmp/simulator/',subj);
dwiFile   = fullfile(rootDir,'dwi_data_b2000_aligned_trilin.nii.gz');
t1File    = fullfile(rootDir,'anatomy','T1w_acpc_dc_restore_1p25.nii.gz');
bvecsFile = fullfile(rootDir,'dwi_data_b2000_aligned_trilin.bvecs');
bvalsFile = fullfile(rootDir,'dwi_data_b2000_aligned_trilin.bvals');
bvecs     = dlmread(bvecsFile);
bvals     = dlmread(bvalsFile);

% load fiber groups
fgArcClean     = fgRead(fullfile(saveDir,'fascicles', strcat('fgrArcClean',subj,'.mat')));
fgCortClean    = fgRead(fullfile(saveDir,'fascicles', strcat('fgrCortClean',subj,'.mat')));
fgArcCortClean = fgRead(fullfile(saveDir,'fascicles', strcat('fgrArcCortClean',subj,'.mat')));

%% load fe structures

feOrig = load(fullfile(saveDir,'fe','fe_rArCort_simulator_105115.mat'));
fe_a99 = load(fullfile(saveDir,'fe','fe_rArcCort_simulator_a0.99_105115.mat'));
fe_a01 = load(fullfile(saveDir,'fe','fe_rArcCort_simulator_a0.01_105115.mat'));
fe_a5  = load(fullfile(saveDir,'fe','fe_rArcCort_simulator_a0.5_105115.mat'));
fe_a999dist = load(fullfile(saveDir,'fe','fe_rArcCort_simulator_a0.999dist_105115.mat'));
fe_a001dist  = load(fullfile(saveDir,'fe','fe_rArcCort_simulator_a0.001dist_105115.mat'));
fe_a5dist   = load(fullfile(saveDir,'fe','fe_rArcCort_simulator_a0.5dist_105115.mat'));

% get new LiFE weights
weights_orig = feGet(feOrig,'fiber weights');
weights_arc  = weights_orig(1:length(fgArcClean.fibers));
weights_cort = weights_orig(1+length(fgArcClean.fibers):end);

new_weights_a99 = feGet(fe_a99,'fiber weights');
new_weights_a01 = feGet(fe_a01, 'fiber weights');
new_weights_a5  = feGet(fe_a5, 'fiber weights');
new_weights_a999dist = feGet(fe_a999dist,'fiber weights');
new_weights_a001dist = feGet(fe_a001dist, 'fiber weights');
new_weights_a5dist   = feGet(fe_a5dist, 'fiber weights');

new_weights_arc999dist = new_weights_a999dist(1:length(weights_arc));
new_weights_arc5dist   = new_weights_a5dist(1:length(weights_arc));
new_weights_arc001dist = new_weights_a001dist(1:length(weights_arc));

new_weights_cort999dist = new_weights_a999dist(1+length(weights_arc):end);
new_weights_cort5dist   = new_weights_a5dist(1+length(weights_arc):end);
new_weights_cort001dist = new_weights_a001dist(1+length(weights_arc):end);


%% get old fiber weights
% set fiber weights using parameters mu, alpha, and beta
mu_arc  = 0.9;  %10^(median(log(weights_arc)));
mu_cort = 0.9;  %10^(median(log(weights_cort)));

% alpha 0.99
alpha99   = 0.99;
beta99    = 1-alpha99;

weight_arc99  = alpha99*mu_arc*ones(size(weights_arc));
weight_cort99 = beta99*mu_cort*ones(size(weights_cort));
weights_a99   = [weight_arc99; weight_cort99];

% alpha 0.5
alpha5  = 0.5;
beta5    = 1-alpha5;

weight_arc5  = alpha5*mu_arc*ones(size(weights_arc));
weight_cort5 = beta5*mu_cort*ones(size(weights_cort));
weights_a5   = [weight_arc5; weight_cort5];

% alpha 0.01
alpha01   = 0.01;
beta01    = 1-alpha01;

weight_arc01  = alpha01*mu_arc*ones(size(weights_arc));
weight_cort01 = beta01*mu_cort*ones(size(weights_cort));
weights_a01   = [weight_arc01; weight_cort01];

% set distributions
% alpha 0.001
alpha   = 0.001;
beta    = 1-alpha;

weight_arc001dist  = alpha*weights_arc;
weight_cort001dist = beta*weights_cort;
weights_a001dist   = [weight_arc001dist; weight_cort001dist];

% alpha 0.5
alpha   = 0.5;
beta    = 1-alpha;

weight_arc5dist  = alpha*weights_arc;
weight_cort5dist = beta*weights_cort;
weights_a5dist   = [weight_arc5dist; weight_cort5dist];

% alpha 0.999
alpha   = 0.999;
beta    = 1-alpha;

weight_arc999dist  = alpha*weights_arc;
weight_cort999dist = beta*weights_cort;
weights_a999dist   = [weight_arc999dist; weight_cort999dist];

%% PLOTS

% original weight distributions 
bins = linspace(-6, 0, 75);
figure('name','Original Weight Distributions','color','w'); 
hold on
[y, x] = hist(log10(weights_arc(weights_arc>0)), bins);
plot(x,y,'-','linewidth',2,'color',[.8 .6 .4])
hold on
[y, x] = hist(log10(weights_cort(weights_cort>0)), bins);   
plot(x, y,'-.','linewidth',2,'color',[.3 .5 .9])
[y, x] = hist(log10(weights_orig(weights_orig>0)), bins);
plot(x, y,'-','linewidth',2,'color',[.3 .8 .4]);

legend({'Arcuate Dist of weights on Y Orig';'Cortico Dist of weights on Y Orig';...
    'Overall Dist of weights'})
set(legend,'Interpreter','latex','fontsize',12)
title('LiFE Weights Test')
ylabel('Number of Fascicles','fontsize',15)
xlabel('log_{10}(Fiber Weight)','fontsize',15)
set(gca,'xlim',[-6 0], 'xtick',[-6 -3 0],'ylim',[ 0 75], 'ytick',...
    [0 37.5 75], 'tickdir', 'out', 'box', 'off', 'fontsize', 15, 'visible', 'on')

% comparison of weight distributions 
bins = linspace(-6, 0, 60);
figure('name','Simulated weights','color','w'); 
hold on
[y, x] = hist(log10(new_weights_a5dist(new_weights_a5dist>0)), bins);
plot(x,y,'-','linewidth',2,'color',[.8 .6 .4])
hold on
[y, x] = hist(log10(new_weights_a001dist(new_weights_a001dist>0)), bins);
plot(x,y,'-','linewidth',2,'color',[.3 .8 .4])

[y, x] = hist(log10(new_weights_a999dist(new_weights_a999dist>0)), bins);   
plot(x, y,'-','linewidth',2,'color',[.3 .5 .9])

[y, x] = hist(log10(weights_a999dist(weights_a999dist>0)), bins);
plot(x, y,'-.','linewidth',2,'color',[.3 .5 .9]);
[y, x] = hist(log10(weights_a001dist(weights_a001dist>0)), bins);
plot(x, y,'-.','linewidth',2,'color',[.3 .8 .4]);
[y, x] = hist(log10(weights_a5dist(weights_a5dist>0)), bins);
plot(x, y,'-.','linewidth',2,'color',[.8 .6 .4]);

legend({'Dist of weights on Y ${\alpha}$ = 0.5';'Dist of weights on Y ${\alpha}$ = 0.001';...
    'Dist of weights on Y ${\alpha}$ = 0.999';'Original weights ${\alpha}$ = 0.999';...
    'Original weights ${\alpha}$ = 0.001';'Original weights ${\alpha}$ = 0.5'})
set(legend,'Interpreter','latex','fontsize',12)
title('Simulated Weights Test')
ylabel('Number of Fascicles','fontsize',15)
xlabel('log_{10}(Fiber Weight)','fontsize',15)
set(gca,'xlim',[-6 0], 'xtick',[-6 -3 0],'ylim',[ 0 75], 'ytick',...
    [0 37.5 75], 'tickdir', 'out', 'box', 'off', 'fontsize', 15, 'visible', 'on')


% all set weights plot
bins = linspace(-6, 0, 60);
figure('name','Set alpha Weights','color','w')
[y,x] = hist(log10(weights_a99(weights_a99>0)),bins);
plot(x,y,'g','linewidth',2)
hold on


[y,x] = hist(log10(weights_a5(weights_a5>0)),bins);
plot(x,y,'k','linewidth',2)

[y,x] = hist(log10(weights_a01(weights_a01>0)),bins);
plot(x,y,'m','linewidth',2)

[y,x] = hist(log10(new_weights_a99(new_weights_a99>0)),bins);
plot(x,y,'b','linewidth',2)

[y,x] = hist(log10(new_weights_a5(new_weights_a5>0)),bins);
plot(x,y,'c','linewidth',2)

[y,x] = hist(log10(new_weights_a01(new_weights_a01>0)),bins);
plot(x,y,'r','linewidth',2)
legend({'Original Set Weights ${\alpha}$=0.99';'Original Set Weights ${\alpha}$=0.5';...
    'Original Set Weights ${\alpha}$=0.01';'LiFE Weights at ${\alpha}$=0.99';...
    'LiFE Weights at ${\alpha}$=0.5';'LiFE Weights at ${\alpha}$=0.01'})
set(legend,'Interpreter','latex','fontsize',12)
title('Weights from Simulated Data and Original Fibers')
ylabel('Number of Fascicles','fontsize',15)
xlabel('log_{10}(Fascicle Weight)','fontsize',15)
set(gca,'xlim',[-6 0], 'xtick',[-6 -3 0],'ylim',[ 0 50], 'ytick',...
    [0 50 100], 'tickdir', 'out', 'box', 'off', 'fontsize', 15, 'visible', 'on')

% all scaled weights plot
bins = linspace(-6, 0, 60);
figure('name','Alpha-scaled Weights','color','w')
[y,x] = hist(log10(weights_a999dist(weights_a999dist>0)),bins);
plot(x,y,'g','linewidth',2)
hold on

[y,x] = hist(log10(weights_a5dist(weights_a5dist>0)),bins);
plot(x,y,'k','linewidth',2)

[y,x] = hist(log10(weights_a001dist(weights_a001dist>0)),bins);
plot(x,y,'m','linewidth',2)

[y,x] = hist(log10(new_weights_a999dist(new_weights_a999dist>0)),bins);
plot(x,y,'b','linewidth',2)

[y,x] = hist(log10(new_weights_a5dist(new_weights_a5dist>0)),bins);
plot(x,y,'c','linewidth',2)

[y,x] = hist(log10(new_weights_a001dist(new_weights_a001dist>0)),bins);
plot(x,y,'r','linewidth',2)
legend({'Original scaled Weights ${\alpha}$=0.99';'Original scaled Weights ${\alpha}$=0.5';...
    'Original scaled Weights ${\alpha}$=0.01';'LiFE Weights at ${\alpha}$=0.99';...
    'LiFE Weights at ${\alpha}$=0.5';'LiFE Weights at ${\alpha}$=0.01'})
set(legend,'Interpreter','latex','fontsize',12)
title('Weights from Simulated Data and Original Fibers')
ylabel('Number of Fascicles','fontsize',15)
xlabel('log_{10}(Fascicle Weight)','fontsize',15)
set(gca,'xlim',[-6 0], 'xtick',[-6 -3 0],'ylim',[ 0 50], 'ytick',...
    [0 50 100], 'tickdir', 'out', 'box', 'off', 'fontsize', 15, 'visible', 'on')



% Arcuate scaled weights plot
bins = linspace(-6, 0, 60);
figure('name','Scaled alpha Arcuate weights','color','w')
[y,x] = hist(log10(weight_arc999dist(weight_arc999dist>0)),bins);
plot(x,y,'b','linewidth',2)
hold on
[y,x] = hist(log10(weight_arc5dist(weight_arc5dist>0)),bins);
plot(x,y,'c','linewidth',2)
[y,x] = hist(log10(weight_arc001dist(weight_arc001dist>0)),bins);
plot(x,y,'r','linewidth',2)
[y,x] = hist(log10(new_weights_arc999dist(new_weights_arc999dist>0)),bins);
plot(x,y,'g','linewidth',2)
[y,x] = hist(log10(new_weights_arc5dist(new_weights_arc5dist>0)),bins);
plot(x,y,'k','linewidth',2)
[y,x] = hist(log10(new_weights_arc001dist(new_weights_arc001dist>0)),bins);
plot(x,y,'m','linewidth',2)
legend({'Original Weights scaled by ${\alpha}$=0.999';'Original Weights scaled by ${\alpha}$=0.5';...
    'Original Weights scaled by ${\alpha}$=0.001';'LiFE Weights at ${\alpha}$=0.999';...
    'LiFE Weights at ${\alpha}$=0.5';'LiFE Weights at ${\alpha}$=0.001'})
set(legend,'Interpreter','latex','fontsize',12)
title('Weights from Simulated Data and Original Fibers')
ylabel('Number of Fascicles','fontsize',15)
xlabel('log_{10}(Fascicle Weight)','fontsize',15)
set(gca,'xlim',[-6 0], 'xtick',[-6 -3 0],'ylim',[ 0 50], 'ytick',...
    [0 50 100], 'tickdir', 'out', 'box', 'off', 'fontsize', 15, 'visible', 'on')


% cortico scaled weights plot
bins = linspace(-6, 0, 60);
figure('name','Scaled alpha Cortico weights','color','w')
[y,x] = hist(weight_cort999dist,bins);
plot(x,y,'b','linewidth',2)
hold on
[y,x] = hist(weight_cort5dist,bins);
plot(x,y,'c','linewidth',2)
[y,x] = hist(weight_cort001dist,bins);
plot(x,y,'r','linewidth',2)
[y,x] = hist(weights_a999dist(1+length(weights_arc):end),bins);
plot(x,y,'g','linewidth',2)
[y,x] = hist(weights_a5dist(1+length(weights_arc):end),bins);
plot(x,y,'k','linewidth',2)
[y,x] = hist(weights_a001dist(1+length(weights_arc):end),bins);
plot(x,y,'m','linewidth',2)
legend({'Original Weights scaled by ${\alpha}$=0.999';'Original Weights scaled by ${\alpha}$=0.5';...
    'Original Weights scaled by ${\alpha}$=0.001';'LiFE Weights at ${\alpha}$=0.999';...
    'LiFE Weights at ${\alpha}$=0.5';'LiFE Weights at ${\alpha}$=0.001'})
set(legend,'Interpreter','latex','fontsize',12)
%axis([ 10^-4.5 10^-.5 -10 1200])
title('Cortico Weights from Simulated Data and Original Fibers')
set(gca,'xscale', 'log')%'Tickdir', 'out', 'Ticklength', [0.01 0.01])

% Arcuate set weights plot
figure('name','Set alpha Arcuate weights','color','w')
[y,x] = hist(weight_arc99,bins);
plot(x,y,'b','linewidth',2)
hold on
[y,x] = hist(weight_arc5,bins);
plot(x,y,'c','linewidth',2)
[y,x] = hist(weight_arc01,bins);
plot(x,y,'r','linewidth',2)
[y,x] = hist(new_weights_a99(1:length(weights_arc)),bins);
plot(x,y,'g','linewidth',2)
[y,x] = hist(new_weights_a5(1:length(weights_arc)),bins);
plot(x,y,'k','linewidth',2)
[y,x] = hist(new_weights_a01(1:length(weights_arc)),bins);
plot(x,y,'m','linewidth',2)
legend({'Original Set Weights ${\alpha}$=0.99';'Original Set Weights ${\alpha}$=0.5';...
    'Original Set Weights ${\alpha}$=0.01';'LiFE Weights at ${\alpha}$=0.99';...
    'LiFE Weights at ${\alpha}$=0.5';'LiFE Weights at ${\alpha}$=0.01'})
set(legend,'Interpreter','latex','fontsize',12)
%axis([ 10^-3 1 -10 1200])
title('Arcuate Weights from Simulated Data and Original Fibers')
set(gca,'xscale', 'log')%'Tickdir', 'out', 'Ticklength', [0.01 0.01])

% Cortico set weights plot
figure('name','Set alpha Cortico weights','color','w')
[y,x] = hist(weight_cort99,bins);
plot(x,y,'b','linewidth',2)
hold on
[y,x] = hist(weight_cort5,bins);
plot(x,y,'c','linewidth',2)
[y,x] = hist(weight_cort01,bins);
plot(x,y,'r','linewidth',2)
[y,x] = hist(new_weights_a99(1+length(weights_arc):end),bins);
plot(x,y,'g','linewidth',2)
[y,x] = hist(new_weights_a5(1+length(weights_arc):end),bins);
plot(x,y,'k','linewidth',2)
[y,x] = hist(new_weights_a01(1+length(weights_arc):end),bins);
plot(x,y,'m','linewidth',2)
legend({'Original Set Weights ${\alpha}$=0.99';'Original Set Weights ${\alpha}$=0.5';...
    'Original Set Weights ${\alpha}$=0.01';'LiFE Weights at ${\alpha}$=0.99';...
    'LiFE Weights at ${\alpha}$=0.5';'LiFE Weights at ${\alpha}$=0.01'})
set(legend,'Interpreter','latex','fontsize',12)
%axis([ 10^-3 1 -10 1200])
title('Cortico Weights from Simulated Data and Original Fibers')
set(gca,'xscale', 'log')%'Tickdir', 'out', 'Ticklength', [0.01 0.01])



% plot
figure
bins = 25;
h1 = histogram(-log(weights_orig),bins,'Normalization','probability');
hold on
h2 = histogram(-log(weights_a999dist),bins,'Normalization','probability');
h3 = histogram(-log(weights_a5dist),bins,'Normalization','probability');
h4 = histogram(-log(weights_a001dist),bins,'Normalization','probability');
h2.FaceColor = 'r';
h3.FaceColor = 'g';
h4.FaceColor = 'c';
h2.BinWidth = h1.BinWidth;
h3.BinWidth = h2.BinWidth;
h4.BinWidth = h3.BinWidth;

figure(5)
plot(h1.BinEdges(1:h1.NumBins)+0.5*h1.BinWidth,h1.Values,'b',...
    h2.BinEdges(1:h2.NumBins)+0.5*h2.BinWidth,h2.Values,'r',...
    h3.BinEdges(1:h3.NumBins)+0.5*h3.BinWidth,h3.Values,'g',...
    h4.BinEdges(1:h4.NumBins)+0.5*h4.BinWidth,h4.Values,'c','linewidth',2)

figure
bins = 25;
h1 = histogram(-log(weights_orig),bins,'Normalization','probability');
hold on
h2 = histogram(-log(new_weights_a99),bins,'Normalization','probability');
h3 = histogram(-log(new_weights_a5),bins,'Normalization','probability');
h4 = histogram(-log(new_weights_a01),bins,'Normalization','probability');
h2.FaceColor = 'r';
h3.FaceColor = 'g';
h4.FaceColor = 'c';
h2.BinWidth = h1.BinWidth;
h3.BinWidth = h2.BinWidth;
h4.BinWidth = h3.BinWidth;

figure
plot(h1.BinEdges(1:h1.NumBins)+0.5*h1.BinWidth,h1.Values,'b',...
    h2.BinEdges(1:h2.NumBins)+0.5*h2.BinWidth,h2.Values,'r',...
    h3.BinEdges(1:h3.NumBins)+0.5*h3.BinWidth,h3.Values,'g',...
    h4.BinEdges(1:h4.NumBins)+0.5*h4.BinWidth,h4.Values,'c','linewidth',2)


%%%%% plots

figure
[y,x] = hist(weights_arc,bins);
plot(x,y,'g','linewidth',2)
hold on
[y,x] = hist(weight_arc001dist,bins);
plot(x,y,'r','linewidth',2)
[y,x] = hist(weight_arc5dist,bins);
plot(x,y,'b','linewidth',2)
[y,x] = hist(weight_arc999dist,bins);
plot(x,y,'c','linewidth',2)


figure
[y,x] = hist(weights_cort,bins);
plot(x,y,'g','linewidth',2)
hold on
[y,x] = hist(weight_cort001dist,bins);
plot(x,y,'r','linewidth',2)
[y,x] = hist(weight_cort5dist,bins);
plot(x,y,'b','linewidth',2)
[y,x] = hist(weight_cort999dist,bins);
plot(x,y,'c','linewidth',2)

figure
[y,x] = hist(weights_arc,bins);
plot(x,y,'g','linewidth',2)
hold on
[y,x] = hist(weight_arc01,bins);
plot(x,y,'r','linewidth',2)
[y,x] = hist(weight_arc5,bins);
plot(x,y,'b','linewidth',2)
[y,x] = hist(weight_arc99,bins);
plot(x,y,'c','linewidth',2)


figure
[y,x] = hist(weights_cort,bins);
plot(x,y,'g','linewidth',2)
hold on
[y,x] = hist(weight_cort01,bins);
plot(x,y,'r','linewidth',2)
[y,x] = hist(weight_cort5,bins);
plot(x,y,'b','linewidth',2)
[y,x] = hist(weight_cort99,bins);
plot(x,y,'c','linewidth',2)

figure
a1 = histogram(weights_arc,bins,'Normalization','probability');
hold on 
a2 = histogram(weight_arc999dist,bins,'Normalization','probability');
a3 = histogram(weight_arc5dist,bins,'Normalization','probability');
a4 = histogram(weight_arc001dist,bins,'Normalization','probability');
a2.BinWidth = a1.BinWidth;
a3.BinWidth = a2.BinWidth;
a4.BinWidth = a3.BinWidth;
a1.FaceColor = 'r';
a2.FaceColor = 'g';
a3.FaceColor = 'c';
a4.FaceColor = 'b';



% find nonzero Arcuate fibers and weights for each alpha
fe_arc_a99 = fe_a99;
zeroArcInd_a99 = find(fe_arc_a99.life.fit.weights(1:length(fgArcClean.fibers))== 0);
nzArcInd_a99   = find(fe_arc_a99.life.fit.weights(1:length(fgArcClean.fibers))~= 0);
fe_arc_a99.fg.fibers = fe_arc_a99.fg.fibers(1:length(fgArcClean.fibers));
fe_arc_a99.fg.fibers(zeroArcInd_a99)=[];
weights_arc_a99 = fe_arc_a99.life.fit.weights(1:length(fgArcClean.fibers));
weights_arc_nz_a99 = weights_arc_a99(nzArcInd_a99);

fe_arc_a5 = fe_a5;
zeroArcInd_a5 = find(fe_arc_a5.life.fit.weights(1:length(fgArcClean.fibers))== 0);
nzArcInd_a5   = find(fe_arc_a5.life.fit.weights(1:length(fgArcClean.fibers))~= 0);
fe_arc_a5.fg.fibers = fe_arc_a5.fg.fibers(1:length(fgArcClean.fibers));
fe_arc_a5.fg.fibers(zeroArcInd_a5)=[];
weights_arc_a5 = fe_arc_a5.life.fit.weights(1:length(fgArcClean.fibers));
weights_arc_nz_a5 = weights_arc_a5(nzArcInd_a5);

fe_arc_a01 = fe_a01;
zeroArcInd_a01 = find(fe_arc_a01.life.fit.weights(1:length(fgArcClean.fibers))== 0);
nzArcInd_a01   = find(fe_arc_a01.life.fit.weights(1:length(fgArcClean.fibers))~= 0);
fe_arc_a01.fg.fibers = fe_arc_a01.fg.fibers(1:length(fgArcClean.fibers));
fe_arc_a01.fg.fibers(zeroArcInd_a01)=[];
weights_arc_a01 = fe_arc_a01.life.fit.weights(1:length(fgArcClean.fibers));
weights_arc_nz_a01 = weights_arc_a01(nzArcInd_a01);

% find nonzero Cortico fibers for each alpha
fe_cort_a99 = fe_a99;
zeroCortInd_a99 = find(fe_cort_a99.life.fit.weights(1+length(fgArcClean.fibers):end)== 0);
nzCortInd_a99 = find(fe_cort_a99.life.fit.weights(1+length(fgArcClean.fibers):end)~= 0);
fe_cort_a99.fg.fibers = fe_cort_a99.fg.fibers(1+length(fgArcClean.fibers):end);
fe_cort_a99.fg.fibers(zeroCortInd_a99)=[];
weights_cort_a99 = fe_cort_a99.life.fit.weights(1+length(fgArcClean.fibers):end);
weights_cort_nz_a99 = weights_cort_a99(nzCortInd_a99);

fe_cort_a5 = fe_a5;
zeroCortInd_a5 = find(fe_cort_a5.life.fit.weights(1+length(fgArcClean.fibers):end)== 0);
nzCortInd_a5 = find(fe_cort_a5.life.fit.weights(1+length(fgArcClean.fibers):end)~= 0);
fe_cort_a5.fg.fibers = fe_cort_a5.fg.fibers(1+length(fgArcClean.fibers):end);
fe_cort_a5.fg.fibers(zeroCortInd_a5)=[];
weights_cort_a5 = fe_cort_a5.life.fit.weights(1+length(fgArcClean.fibers):end);
weights_cort_nz_a5 = weights_cort_a5(nzCortInd_a5);

fe_cort_a01 = fe_a01;
zeroCortInd_a01 = find(fe_cort_a01.life.fit.weights(1+length(fgArcClean.fibers):end)== 0);
nzCortInd_a01 = find(fe_cort_a01.life.fit.weights(1+length(fgArcClean.fibers):end)~= 0);
fe_cort_a01.fg.fibers = fe_cort_a01.fg.fibers(1+length(fgArcClean.fibers):end);
fe_cort_a01.fg.fibers(zeroCortInd_a01)=[];
weights_cort_a01 = fe_cort_a01.life.fit.weights(1+length(fgArcClean.fibers):end);
weights_cort_nz_a01 = weights_cort_a01(nzCortInd_a01);

weights_nz_a01 = [weights_arc_nz_a01; weights_cort_nz_a01];
weights_nz_a5  = [weights_arc_nz_a5 ; weights_cort_nz_a5 ];
weights_nz_a99 = [weights_arc_nz_a99; weights_cort_nz_a99];

%% The numbers
% alpha 0.99: 107 nonzero out of 760 Cortico fibers, 66 nonzero out of 470
% Arcuate fibers
avg_arc_a99  = mean(weights_arc_nz_a99);  % 0.0354
avg_cort_a99 = mean(weights_cort_nz_a99); % 0.0388

% alpha 0.5: 96 nonzero out of 760 Cortico fibers, 67 nonzero out of 470
% Arcuate fibers
avg_arc_a5  = mean(weights_arc_nz_a5);   % 0.0279
avg_cort_a5 = mean(weights_cort_nz_a5);  % 0.0328

% alpha 0.01: 104 out of 760 Cortico fibers, 65 nonzero out of 470 Arcuate
% fibers 
avg_arc_a01  = mean(weights_arc_nz_a01);  % 0.0297
avg_cort_a01 = mean(weights_cort_nz_a01); % 0.0335

% remove zero-weighted fibers
zeroInd_a99 = find(new_weights_a99==0);
fe_a99_nz = fe_a99;
fe_a99_nz.fg.fibers(zeroInd_a99)=[];

zeroInd_a5 = find(new_weights_a5==0);
fe_a5_nz = fe_a5;
fe_a5_nz.fg.fibers(zeroInd_a5)=[];

zeroInd_a01 = find(new_weights_a01==0);
fe_a01_nz = fe_a01;
fe_a01_nz.fg.fibers(zeroInd_a01)=[];

% plot fiber distributions

%weight_bins = log10(logspace(10^-6,10^-5,60));
weight_bins = log10(logspace(min(weights_arc_nz_a99),max(weights_arc_nz_a99)));

figure (1)
%subplot(1,2,1)

[y,x] = hist(weights_a99(1:length(weights_arc)),weight_bins);
plot(x,y,'-r')


hold on
[y,x] = hist(weights_arc_nz_a99,weight_bins); 
plot(x,y,'-b')
legend('Original, set weights for Arcuate fibers alpha 0.99',...
    'Nonzero LiFE weights for Arcuate alpha 0.99','location','northeast')
axis([0 0.3 0 10]);

%subplot(1,2,2)
hist(weights_a99(1+length(weights_arc):end))
title('Original, set weights for Cortico fibers beta 0.01')

subplot(2,2,4)
hist(weights_cort_nz_a99); 
title('Nonzero LiFE weights for Cortico fibers beta 0.01')


figure (2)
subplot(2,2,1)
hist(weights_a5(1:length(weights_arc)))
title('Original, set weights for Arcuate fibers alpha 0.5')

subplot(2,2,2)
hist(weights_arc_nz_a5)
title('Nonzero LiFE weights for Arcuate fibers alpha 0.5')

subplot(2,2,3)
hist(weights_a5(1+length(weights_arc):end))
title('Original, set weights for Cortico fibers beta 0.5')

subplot(2,2,4)
hist(weights_cort_nz_a5)
title('Nonzero LiFE weights for Cortico fibers beta 0.5')

figure (3)
subplot(2,2,1)
hist(weights_a01(1:length(weights_arc)))
title('Original, set weights for Arcuate fibers alpha 0.01')

subplot(2,2,2)
hist(weights_arc_nz_a01)
title('Nonzero LiFE weights for Arcuate fibers alpha 0.01')

subplot(2,2,3)
hist(weights_a01(1+length(weights_arc):end))
title('Original, set weights for Cortico fibers beta 0.99')

subplot(2,2,4)
hist(weights_cort_nz_a01)
title('Nonzero LiFE weights for Cortico fibers beta 0.99')

figure (4)
subplot(1,3,1)
for ii = 1:length(fe_arc_a99.fg.fibers);
    plot3(fe_arc_a99.fg.fibers{ii}(1,:),fe_arc_a99.fg.fibers{ii}(2,:),fe_arc_a99.fg.fibers{ii}(3,:),'b'); hold on
end
hold on
for ii = 1:length(fe_cort_a99.fg.fibers);
    plot3(fe_cort_a99.fg.fibers{ii}(1,:),fe_cort_a99.fg.fibers{ii}(2,:),fe_cort_a99.fg.fibers{ii}(3,:),'r'); hold on
end
view(90,0)
title('Nonzero Arcuate and Cortico Fibers alpha 0.99')



subplot(1,3,2)
for ii = 1:length(fe_arc_a5.fg.fibers);
    plot3(fe_arc_a5.fg.fibers{ii}(1,:),fe_arc_a5.fg.fibers{ii}(2,:),fe_arc_a5.fg.fibers{ii}(3,:),'b'); hold on
end
hold on
for ii = 1:length(fe_cort_a5.fg.fibers);
    plot3(fe_cort_a5.fg.fibers{ii}(1,:),fe_cort_a5.fg.fibers{ii}(2,:),fe_cort_a5.fg.fibers{ii}(3,:),'r'); hold on
end
view(90,0)
title('Nonzero Arcuate and Cortico Fibers alpha 0.5')



subplot(1,3,3)
for ii = 1:length(fe_arc_a01.fg.fibers);
    plot3(fe_arc_a01.fg.fibers{ii}(1,:),fe_arc_a01.fg.fibers{ii}(2,:),fe_arc_a01.fg.fibers{ii}(3,:),'b'); hold on
end
hold on
for ii = 1:length(fe_cort_a01.fg.fibers);
    plot3(fe_cort_a01.fg.fibers{ii}(1,:),fe_cort_a01.fg.fibers{ii}(2,:),fe_cort_a01.fg.fibers{ii}(3,:),'r'); hold on
end
view(90,0)
title('Nonzero Arcuate and Cortico Fibers alpha 0.01')


bins = 25;
figure('name','Weights Distribution','color','w');
h = histogram(-log(orig_weights),bins,'Normalization','probability');
hold on
h1 = histogram(-log(orig_arc_weights),bins,'Normalization','probability');
h2 = histogram(-log(orig_cort_weights),bins,'Normalization','probability');

h1.BinWidth = h.BinWidth;
h2.BinWidth = h1.BinWidth;

h.FaceColor = 'r';
h1.FaceColor = 'b';
h2.FaceColor = 'g';


plot(h1.BinEdges(1:h1.NumBins)+0.5*h1.BinWidth,h1.Values,'r',...
    h2.BinEdges(1:h2.NumBins)+0.5*h2.BinWidth,h2.Values,'g',...
    h.BinEdges(1:h.NumBins)+0.5*h.BinWidth,h.Values,'b','linewidth',2)

