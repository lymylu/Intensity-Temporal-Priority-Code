%%
% THE tensormatrix_initial is the matrixs of channel*time*freq*(trial*subject), rank is the number of components
%tensormatrix_initial=(tensormatrix_initial-repmat(min(reshape(tensormatrix_initial,[],1,size(tensormatrix_initial,3)),[],1),[size(tensormatrix_initial,1),size(tensormatrix_initial,2),1]))./...
%    (repmat((max(reshape(tensormatrix_initial,[],1,size(tensormatrix_initial,3)),[],1)-min(reshape(tensormatrix_initial,[],1,size(tensormatrix_initial,3)),[],1)),[size(tensormatrix_initial,1),size(tensormatrix_initial,2),1]));
% % tensormatrix were normalized across each neuron.
%% transfer the tensor1.h5 result from TCA_cal.py to ktensor
clear all;
preprocess={'_nose','_Fz','_averageref'};
all=2;
basepath='/mnt/Share/yuelp/LEP_4L_3s_IClabelmoreremove/Intensitypriority';
tensormat=matfile(fullfile(basepath,['trialaverage',preprocess{all},'.mat']),'Writable',true);
tensormatrix_cal=tensormat.tensormatrix_mwt; 
t=linspace(-1,1.999,3000);

tensormatrix_cal=tensor(permute(tensormatrix_cal(t>=-0.2&t<=0.5,:,:,:),[3,1,2,4]));
t=t(t>=-0.2&t<=0.5);
Datasets=h5info(fullfile(basepath,['tensor_mwt',preprocess{all},'_silice.h5'])); % tensor_mwt_averageref or tensor_freqaverage_mwt
tensorfilename=fullfile(basepath,['tensor_mwt',preprocess{all},'_silice.h5']);
tensor_py=[];rmse_py=[];similarity_py=[];
rank=30;
for i=1:rank
    lambda=h5read(tensorfilename,['/factor_lam_',num2str(i)]);
    factorchannel=h5read(tensorfilename,['/factor_channel_',num2str(i)]);
    factortime=h5read(tensorfilename,['/factor_time_',num2str(i)]);
    factortrial=h5read(tensorfilename,['/factor_trial_',num2str(i)]);
    factorfreq=h5read(tensorfilename,['/factor_freq_',num2str(i)]);
    rmse_tmp=h5read(tensorfilename,['/reconerror_',num2str(i)]);
    rmse_py(i,:)=rmse_tmp;
    similarity_tmp=h5read(tensorfilename,['/similarity_',num2str(i)]);
    similarity_py(i,:)=similarity_tmp;
    for c=1:10
        tensor_py{i,c}=ktensor(squeeze(lambda(:,c)),squeeze(factorchannel(:,:,c))',squeeze(factortime(:,:,c))',squeeze(factorfreq(:,:,c))',squeeze(factortrial(:,:,c))');
    end
end
% calculate similarity and rmse
for R=1:rank
    for iter=1:10
        rmse(R,iter)=fg(tensor_py{R,iter},tensormatrix_cal,norm(tensormatrix_cal).^2,norm(tensormatrix_cal).^2,false,false,false);
        similarity(R,iter)=score(tensor_py{R,iter},tensor_py{R,1});
    end
end
figure; subplot(2,2,1);scatter(1:rank,rmse_py); hold on; plot(1:rank,mean(rmse_py,2)); title('reconstructed error of different TCA model rank');
subplot(2,2,2);scatter(2:rank,diff(rmse_py)); hold on; plot(2:rank,mean(diff(rmse_py),2)); title('Difference of the reconstructed error');
subplot(2,2,3);scatter(1:rank,similarity_py(:,2:end)); hold on; plot(1:rank,mean(similarity_py(:,2:10),2)); title('similarity of different TCA model rank');
% svm decode accuarray using Tensor_all to identify the best rank
eventlistall=cat(1,repmat({'Righthigh'},[95,1]),repmat({'Rightlow'},[95,1]),repmat({'Lefthigh'},[95,1]),repmat({'Leftlow'},[95,1]));
eventnumel=[];
for i=1:4
    eventnumel=cat(1,eventnumel,i*ones(95,1));
end 
for R=1:rank
   for iter=1:10
       mdl=fitcecoc(tensor_py{R,iter}.U{4},eventlistall,'Coding','onevsall');
       label=predict(mdl,tensor_py{R,iter}.U{4});
       accuarray(R,iter)=sum(cellfun(@(x,y) strcmp(x,y),label,eventlistall,'UniformOutput',1))/380;
        %accuarray(R,iter,i)=svmtrain(cate,Tensor_all{R,iter}.U{3},' -s 0 -t 0 -v 5 -h 0');
    end
end
subplot(2,2,4); hold off;scatter(1:rank,mean(accuarray,2)); title('SVM predict performance of all 4 conditions');
% for i=1:10
% [res_x, idx_of_result(i)] = knee_pt(accuarray(:,i)',1:30);
% end
for i=1:10
v=locfit((1:rank)',accuarray(:,i));
accuarraynew(:,i)=predict(v,(1:30)');
end
hold on; plot(1:rank,mean(accuarraynew,2));
idx_of_result=knee_pt(mean(accuarraynew,2)); 

%savefig(gcf,fullfile(basepath,['TCAevaluate',preprocess{all},'.fig']));
%%
%  generate the sampled dataspike matrix (trials were subsampled)
tensormatrix=tensormat.tensormatrix_mwt;
datadeconv_permute=[];index_all=[];
for i=1:20
    group1index=randperm(95,80);
    datadeconv_permute(:,:,:,1:80,i)=tensormatrix(:,:,:,group1index);
    group1index=randperm(95,80);
    datadeconv_permute(:,:,:,81:160,i)=tensormatrix(:,:,:,group1index+95);
     group1index=randperm(95,80);
    datadeconv_permute(:,:,:,161:240,i)=tensormatrix(:,:,:,group1index+190);
    group1index=randperm(95,80);
     datadeconv_permute(:,:,:,241:320,i)=tensormatrix(:,:,:,group1index+285);
end
tensormat.tensormatrix_resample=datadeconv_permute;



%% transfer the tensorpy_resample.h5 to ktensor
Datasets=h5info(fullfile(basepath,['tensor_mwt',preprocess{all},'_resample_13.h5']));
tensorfilename=fullfile(basepath,['tensor_mwt',preprocess{all},'_resample_13.h5']);
tensor_py_resample=[];
for i=1:20
     lambda=h5read(tensorfilename,['/factor_lam_',num2str(i-1)]);
    factorchannel=h5read(tensorfilename,['/factor_channel_',num2str(i-1)]);
    factortime=h5read(tensorfilename,['/factor_time_',num2str(i-1)]);
    factortrial=h5read(tensorfilename,['/factor_trial_',num2str(i-1)]);
    factorfreq=h5read(tensorfilename,['/factor_freq_',num2str(i-1)]);
    for c=1:10
        tensor_py_resample{i,c}=ktensor(squeeze(lambda(:,c)),squeeze(factorchannel(:,:,c))',squeeze(factortime(:,:,c))',squeeze(factorfreq(:,:,c))',squeeze(factortrial(:,:,c))');
    end
end

%% find the rank for TCAdiscetized deconvolved neural activity 
% [Tensor_all,Tensor_origin,similarity,rmse]=getTCArank(tensormatrix_cal);
% save('TCAworkspace.mat','tensormatrix_cal','Tensor_all','Tensor_origin','similarity','rmse','-v7.3');
%%
figure; 
subplot(1,2,1);scatter(1:40,rmse);
subplot(1,2,2);scatter(1:40,similarity(:,2:10));hold on; plot(1:40,mean(similarity(:,2:10),2));
%%
vizopts = {'PlotCommands',{'bar','line','bar','scatter'},...
'ModeTitles',{'Channel','Time','freq','Trialtype'},...
'BottomSpace',0.10,'HorzSpace',0.04,'YLim',{[0,0.2],[0,0.1],[0,1],[0,0.2]},'YTicks',false,'HorzSpace',0.02};
tmp=normalize(tensor_py{12,1});
info1 = viz(normalize(tensor_py{12,1}),'Figure',1,vizopts{:});
for i=1:length(info1.FactorAxes)
     axes(info1.FactorAxes(1,i));
     topoplot(tmp.U{1}(:,i),chanlocs,'emarkersize',4); caxis([0,0.3]);
     %replot_NeuronType(gca,neurontypelabel);
     axes(info1.FactorAxes(4,i));
     replot_TrialType(gca,eventlistall);
end
% figure;a=1;
% for i=1:12
%     subplot(4,5,a)
%      topoplot(tmp.U{1}(:,i),chanlocs,'emarkersize',4); caxis([0,0.3]);
%      a=a+1;
% end
%% end
figure; 
for i=1:length(eventnumber)
    plot(1:40,squeeze(mean(accuarray(:,:,i),2)));
    hold on;
end
%% plot the temperal factor which encode side and intensity.
M1=permute(tensor_py{12,1},[1,2,3,4]); % 12 for nose short, 13 for Fz short % 19 for Fz, 21 for nose
M1=normalize(M1);
%M1.U{3}=basecorrect(M1.U{3},1:6,1,6,'relativepower');
t=linspace(-1,1.999,3000);
t=t(t>=-0.2&t<=0.5);
chanlocs=tensormat.chanlocs;
Intensityinfo=cat(1,repmat({'High'},[190,1]),repmat({'Low'},[190,1]));
Intensityinfo=cat(1,repmat(1,[190,1]),repmat(0,[190,1]));
Sideinfo=repmat(cat(1,repmat({'Right'},[95,1]),repmat({'Left'},[95,1])),[2,1]);
Sideinfo=repmat(cat(1,repmat(1,[95,1]),repmat(0,[95,1])),[2,1]);
[Sidecomponent,Intensitycomponent,Bothcomponent,statp]=TCAfeatureselect(M1,Sideinfo,Intensityinfo);
figure;
subplot(1,4,1)
hold on;
%M1.U{2}=basecorrect(M1.U{2},t,-0.2,0.5,'normalized');
time_intensity=M1.U{2}(:,Intensitycomponent);
time_intensity=basecorrect(mean(time_intensity,2),t',-0.2,0,'subtract');
plot(t,mean(time_intensity,2)); 
time_both=M1.U{2}(:,Bothcomponent);
time_both=basecorrect(mean(time_both,2),t,-0.2,0,'subtract');
plot(t,mean(time_both,2)); axis tight 
xlim([-0.2,0.5]);
try
    time_side=M1.U{2}(:,Sidecomponent);
    time_both=basecorrect(mean(time_both,2),t,-0.2,0,'subtract');
plot(t,mean(time_side,2)); axis tight
end
title('Temporal Weights'); legend({'Intensity','Both factor'});
%M1.U{1}=basecorrect(M1.U{1},1:60,1,60,'normalized');
subplot(1,4,2); topo_both=M1.U{1}(:,Bothcomponent);topoplot(mean(topo_both,2),chanlocs); caxis([0,0.2]); title("Both Topograph")
subplot(1,4,3);topo_intensity=M1.U{1}(:,Intensitycomponent);topoplot(mean(topo_intensity,2),chanlocs); caxis([0,0.2]); title("Intensity Topograph");
subplot(1,4,4);
hold on;
freq=M1.U{3};
%freq=basecorrect(freq,1:6,1,6,'normalized');
freq_Both=freq(:,Bothcomponent);
%freq_side=basecorrect(mean(freq_side,2),1:6,1,6,'normalized');
bar(1:6,mean(freq_Both,2),0.3,'r');
freq_intensity=freq(:,Intensitycomponent);
%freq_intensity=basecorrect(mean(freq_intensity,2),1:6,1,6,'normalized');
bar((1:6)+0.5,mean(freq_intensity,2),0.3,'b'); title('Freq Weights'); legend({'Both factor','Intensity factor'});
savefig(gcf,fullfile(basepath,['TCAresult',preprocess{all},'.fig']));
%%

figure; 
for i=1:9
    channelweights=tensor_py{11,1}.U{1}(:,i);
    subplot(5,5,i); topoplot(channelweights,EEG.chanlocs); caxis([-5,5]);
    if ismember(i,Bothcomponent)
        title('Both');
    elseif ismember(i,Intensitycomponent)
        title('Intensity');
    end
end
%% feature stability across rank
rankrange=[10,11,12,13,14];
bothpeak_rank=[];
intenpeak_rank=[];
freq_inten_rank=[];
freq_both_rank=[];
for r=1:length(rankrange)
M1=permute(tensor_py{rankrange(r),1},[1,2,3,4]); % 12 for nose short, 13 for Fz short % 19 for Fz, 21 for nose
M1=normalize(M1);
%M1.U{3}=basecorrect(M1.U{3},1:6,1,6,'relativepower');
t=linspace(-1,1.999,3000);
t=t(t>=-0.2&t<=0.5);
chanlocs=tensormat.chanlocs;
Intensityinfo=cat(1,repmat({'High'},[190,1]),repmat({'Low'},[190,1]));
Intensityinfo=cat(1,repmat(1,[190,1]),repmat(0,[190,1]));
Sideinfo=repmat(cat(1,repmat({'Right'},[95,1]),repmat({'Left'},[95,1])),[2,1]);
Sideinfo=repmat(cat(1,repmat(1,[95,1]),repmat(0,[95,1])),[2,1]);
[Sidecomponent,Intensitycomponent,Bothcomponent,statp]=TCAfeatureselect(M1,Sideinfo,Intensityinfo);
time_intensity=M1.U{2}(:,Intensitycomponent); 
time_both=M1.U{2}(:,Bothcomponent)
tmppeak=findpeaks(mean(time_intensity,2));
tmp=min(find(tmppeak.loc>60));
intenpeak_rank(r)=t(tmppeak.loc(tmp));
tmppeak=findpeaks(mean(time_both,2));
tmp=min(find(tmppeak.loc>60));
bothpeak_rank(r)=t(tmppeak.loc(tmp));
freq=M1.U{3};
freq_both_rank(:,r)=mean(freq(:,Bothcomponent),2);
freq_inten_rank(:,r)=mean(freq(:,Intensitycomponent),2);
end
figure;plot(rankrange, intenpeak_rank); hold on; plot(rankrange,bothpeak_rank);
figure;surf(freq_inten_rank,repmat(0.5,[size(freq_inten_rank),3]));
hold on;surf(freq_both_rank,repmat(0.8,[size(freq_inten_rank),3]));
%% summarized the factors from tensorpy_resample.h5 and plot the permutated temporal factors which encode intensity or side
t=linspace(-0.2,0.5,701);    
neuronsidecomptmp=[];tempsidecomptmp=[];
neuronintencomptmp=[];tempintencomptmp=[];
Intensityinfo=cat(1,repmat({'High'},[160,1]),repmat({'Low'},[160,1]));
Intensityinfo=cat(1,repmat(1,[160,1]),repmat(0,[160,1]));
Sideinfo=repmat(cat(1,repmat({'Right'},[80,1]),repmat({'Left'},[80,1])),[2,1]);
Sideinfo=repmat(cat(1,repmat(1,[80,1]),repmat(0,[80,1])),[2,1]);
timeintencomptmp=[];timebothcomptmp=[];freqbothcomptmp=[];freqintensitycomptmp=[];
for i=1:20 % resample number =20
    M1=permute(tensor_py_resample{i,1},[1,2,3,4]);
    M1=normalize(M1);
    [Sidecomponent,Intensitycomponent,Bothcomponent]=TCAfeatureselect(M1,Sideinfo,Intensityinfo); 
    time_both=M1.U{2}(:,Bothcomponent);
    %time_both=basecorrect(time_both,t,-0.2,0,'subtract');
     time_intensity=M1.U{2}(:,Intensitycomponent);
   %time_intensity=basecorrect(time_intensity,t,-0.2,0,'subtract');
    timeintencomptmp(:,i)=mean(time_intensity,2);
    timebothcomptmp(:,i)=mean(time_both,2);
    freq_permute=M1.U{3};
     %freq_permute=basecorrect(freq_permute,1:6,1,6,'relativepower');
    freq_both=freq_permute(:,Bothcomponent);
    %freq_both=basecorrect(mean(freq_both,2),1:6,1,6,'relativepower');
    freqbothcomptmp(:,i)=mean(freq_both,2);
    freq_intensity=freq_permute(:,Intensitycomponent);
    %freq_intensity=basecorrect(mean(freq_intensity,2),1:6,1,6,'relativepower');
    freqintensitycomptmp(:,i)=mean(freq_intensity,2);
end
figure; subplot(1,2,1);
timeintencomptmp=basecorrect(timeintencomptmp,t,-0.2,0,'subtract');
timebothcomptmp=basecorrect(timebothcomptmp,t,-0.2,0,'subtract');
plot(t',mean(timeintencomptmp,2),'b');hold on;
plot(t',mean(timebothcomptmp,2),'r');hold on;
shadebar(t',timeintencomptmp,'b');hold on;
shadebar(t',timebothcomptmp,'r');
for i=1:20
    tmppeak=findpeaks(timeintencomptmp(:,i));
   tmp=min(find(tmppeak.loc>60));
    intenpeak(:,i)=t(tmppeak.loc(tmp));
    tmppeak=findpeaks(timebothcomptmp(:,i));
   tmp=min(find(tmppeak.loc>60));
    bothpeak(:,i)=t(tmppeak.loc(tmp));
end
subplot(1,2,2)
errorbar(1:6,mean(freqintensitycomptmp,2),std(freqintensitycomptmp')','b');hold on;
errorbar(1:6,mean(freqbothcomptmp,2),std(freqintensitycomptmp')','r');hold on;
savefig(gcf,fullfile(basepath,'Figure',['TCA_shuffle',preprocess{all},'.fig']));
for i=1:6
    [~,p(i)]=ttest(freqintensitycomptmp(i,:),freqbothcomptmp(i,:));
end
save(fullfile(basepath,'Figure',['TCA_shuffle_stat',preprocess{all},'.mat']),'freqintensitycomptmp','freqbothcomptmp','intenpeak','bothpeak');

function [Sidecomponent,Intensitycomponent,Bothcomponent,statp]=TCAfeatureselect(M1,Sideinfo,Intensityinfo)
% return the component index which significantly related to Side or Intensity.
    trialcomponent=M1.U{4}; % using trial factor to fit.
    trialstat=trialcomponent;
    indexvalid=[];statp=[];
for i=1:size(trialcomponent,2)
    [~,~,stat]=glmfit([Sideinfo,Intensityinfo],trialcomponent(:,i),"normal","constant","on");
    statp(:,i)=stat.p(2:3);
end
statp=statp';
statp=mafdr(statp(:),'BHFDR',true);
statp=reshape(statp,size(M1.U{4},2),[]);
indexvalid=[];
for i=1:2
%     if i==1
%         indexvalid(:,i)=statp(:,1)<0.05&statp(:,2)>0.05;
%     elseif i==2
%         indexvalid(:,i)=statp(:,2)<0.05&statp(:,1)>0.05;
%     elseif i==3
%         indexvalid(:,i)=statp(:,1)<0.05&statp(:,2)<0.05;
%     end
   indexvalid(:,i)=statp(:,i)<0.05;
end
Sidecomponent=find(indexvalid(:,1)==1&indexvalid(:,2)~=1);
Intensitycomponent=find(indexvalid(:,2)==1&indexvalid(:,1)~=1);
Bothcomponent=find(indexvalid(:,2)==1&indexvalid(:,1)==1);
end
function replot_TrialType(figaxes,trialtypelabel)
% replot the current figure to differentiate trial types from the trial
% associate the scatter plot of TCs and reconstructed TCs
 axes(figaxes);
 data=findobj(figaxes);
 Y_data=data(end).YData;
 X_data=data(end).XData;
 delete(data(end));
 trialtypes=unique(trialtypelabel,'stable');
for i=1:length(trialtypes)
scatter(X_data(ismember(trialtypelabel,trialtypes{i})),Y_data(ismember(trialtypelabel,trialtypes{i})),10,'filled')
end
end
function replot_NeuronType(figaxes,celltypelabel)
% replot the current figure to differentiate cell types from the neuron
% associate the bar plot of TCs and reconstructed TCs
 axes(figaxes);
 data=findobj(figaxes);
 Y_data=data(end).YData;
 X_data=data(end).XData;
 delete(data(end));hold on;
 bar(X_data(ismember(celltypelabel,'all')),Y_data(ismember(celltypelabel,'all')),'b','BarWidth',0.5,'EdgeAlpha',0);
 bar(X_data(ismember(celltypelabel,'pv')),Y_data(ismember(celltypelabel,'pv')),'r','BarWidth',0.5,'EdgeAlpha',0);
 
end