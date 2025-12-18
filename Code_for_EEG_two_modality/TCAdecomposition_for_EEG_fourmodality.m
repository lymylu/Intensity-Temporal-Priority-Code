% load tensormatrix and tensor results.
clear all;
tensormat=matfile('/mnt/Share/yuelp/Data_Fourmodalities/Intensitypriority/tensormatrix_mwt.mat');
%tensormat=matfile('/mnt/Share/yuelp/Data_Fourmodalities/Intensitypriority/trialaverage_nose.mat');
% channel* time *frequency * [condition (LH,LL,RH,RL)*subject]*modality
t=linspace(-0.2,0.5,701);
tensormatrix_cal=tensormat.Tensormatrix;
%tensormatrix_origin=tensormat.Specsubject_mwt;
%
% tensormatrix_cal=[];
% for i=1:4
%     tensormatrix_cal=cat(5,tensormatrix_cal,squeeze(tensormatrix_origin(:,:,:,2,i,:)));
% end
% tensormatrix_cal=tensor(tensormatrix_cal);
tensormatrix_cal=tensor(squeeze(tensormatrix_cal(:,:,:,:,1)));
tensorfilename=['/mnt/Share/yuelp/Data_Fourmodalities/Intensitypriority/tensor_laser_mwt.h5'];
%tensorfilename=['/mnt/Share/yuelp/Data_Fourmodalities/Intensitypriority/tensor_elec_only.h5'];

tensor_py=[];rmse_py=[];similarity_py=[];
subjectnum=107;
% type order RightHigh,RightLow,LeftHigh,LeftLow for eeg four modalities.
trialtypelabel=repmat({'Lefthigh','Leftlow','Righthigh','Rightlow'},[subjectnum,1]);
%trialtypelabel=repmat({'Lefthigh','Righthigh','Leftlow','Rightlow'},[subjectnum,1]);
trialtypelabel=reshape(trialtypelabel,[],1);
for i=1:30
    lambda=h5read(tensorfilename,['/factor_lam_',num2str(i)]);
    factortime=h5read(tensorfilename,['/factor_time_',num2str(i)]);
    factortrial=h5read(tensorfilename,['/factor_trial_',num2str(i)]);
    factorchannel=h5read(tensorfilename,['/factor_channel_',num2str(i)]);
    factorfrequency=h5read(tensorfilename,['/factor_freq_',num2str(i)]);
    %factormodality=h5read(tensorfilename,['/factor_modality_',num2str(i)]);
    rmse_tmp=h5read(tensorfilename,['/reconerror_',num2str(i)]);
    rmse_py(i,:)=rmse_tmp;
    similarity_tmp=h5read(tensorfilename,['/similarity_',num2str(i)]);
    similarity_py(i,:)=similarity_tmp;
    for c=1:10
       % tensor_py{i,c}=ktensor(squeeze(lambda(:,c)),squeeze(factorchannel(:,:,c))',squeeze(factortime(:,:,c))',squeeze(factorfrequency(:,:,c))',squeeze(factortrial(:,:,c))',squeeze(factormodality(:,:,c))');
    tensor_py{i,c}=ktensor(squeeze(lambda(:,c)),squeeze(factorchannel(:,:,c))',squeeze(factortime(:,:,c))',squeeze(factorfrequency(:,:,c))',squeeze(factortrial(:,:,c))');
    end
end
figure; subplot(2,2,1);scatter(1:30,rmse_py); hold on; plot(1:30,mean(rmse_py,2)); title('reconstructed error of different TCA model rank'); xlim([2,30]);
subplot(2,2,2);scatter(2:30,diff(rmse_py)); hold on; plot(2:30,mean(diff(rmse_py),2)); title('Difference of the reconstructed error');xlim([2,30]);
subplot(2,2,3);scatter(1:30,similarity_py); hold on; plot(1:30,mean(similarity_py(:,2:10),2)); title('similarity of different TCA model rank');xlim([2,30]);
% svm decode accuarray using Tensor_all to identify the best rank
eventlistall=trialtypelabel;
for R=1:30
   for iter=1:10
       mdl=fitcecoc(tensor_py{R,iter}.U{4},eventlistall,'Coding','onevsall');
       label=predict(mdl,tensor_py{R,iter}.U{4}); 
       accuarray(R,iter)=sum(cellfun(@(x,y) strcmp(x,y),label,eventlistall,'UniformOutput',1))/length(eventlistall);
    end
end
subplot(2,2,4); hold off;scatter(1:30,accuarray(:,1)); title('SVM predict performance of all 4 conditions'); xlim([2,30]);
for i=1:10
v=locfit((1:30)',accuarray(:,i));
accuarraynew(:,i)=predict(v,(1:30)');
end
hold on; plot(1:30,accuarraynew(:,1));
%idx_of_result=knee_pt(accuarraynew(1:end,1)); 
%%
%  generate the sampled dataspike matrix (trials were resampled)
tensormat=matfile('/mnt/Share/yuelp/Data_Fourmodalities/Intensitypriority/tensormatrix_mwt.mat','Writable',true);
tensormatrix=tensormat.Tensormatrix; % channel*time*freq*trial*modality
data_resample=[];index_all=[];
for i=1:20
    data_resample(:,:,:,:,:,i)=cat(4,tensormatrix(:,:,:,datasample(1:107,90,'Replace',false),:),tensormatrix(:,:,:,datasample(108:214,90,'Replace',false),:),tensormatrix(:,:,:,datasample(215:321,90,'Replace',false),:),tensormatrix(:,:,:,datasample(322:428,90,'Replace',false),:));
end
tensormat.Tensormatrix_resample=data_resample;
%% transfer the tensorpy_resample.h5 to ktensor
tensorfilename=['/mnt/Share/yuelp/TCAdecomposition/tensor_eeg_lasermodalities_resample_new.h5'];
tensor_py_resample=[];
for i=1:20
    lambda=h5read(tensorfilename,['/factor_lam_',num2str(i-1)]);
    factortime=h5read(tensorfilename,['/factor_time_',num2str(i-1)]);
    factortrial=h5read(tensorfilename,['/factor_trial_',num2str(i-1)]);
    factorchannel=h5read(tensorfilename,['/factor_channel_',num2str(i-1)]);
    factorfrequency=h5read(tensorfilename,['/factor_freq_',num2str(i-1)]);
    for c=1:10
        tensor_py_resample{i,c}=ktensor(squeeze(lambda(:,c)),squeeze(factorchannel(:,:,c))',squeeze(factortime(:,:,c))',squeeze(factorfrequency(:,:,c))',squeeze(factortrial(:,:,c))');
    end
 end
%%
vizopts = {'PlotCommands',{'bar','line','bar','scatter','bar'},...
'ModeTitles',{'channel','Time','freqeuncy','Trialtype','modality'},...
'BottomSpace',0.10,'HorzSpace',0.04,'Normalize',1,'YLim',{'same','same','same','same','same'},'YTicks',false,'HorzSpace',0.02};
info1 = viz(tensor_py{17,1},'Figure',1,vizopts{:});
M1=tensor_py{17,1};
M1=normalize(M1,1);
figure;
for i=1:17
    subplot(5,5,i);
    topoplot(M1.U{1}(:,i),chanlocs); %caxis([0,0.2]); 
end
%% plot the Laser/Elec induced Time TCs. Both Rank=17
% 1. select the modality specific TCs 
% 2. plot the temperal factor which encode side and intensity.
load('/mnt/Share/yuelp/Data_Fourmodalities/Intensitypriority/trialaverage_nose.mat','chanlocs');
t=linspace(-0.2,0.5,701);

M1=tensor_py{17,1};%t=linspace(-0.1,0.3,400);
M1=normalize(M1);
timepy=M1.U{2}; validindex=[];peakindex=[];

Sideinfo=cat(1,repmat(0,[2*subjectnum,1]),repmat(1,[2*subjectnum,1]));
Intensityinfo=repmat(cat(1,repmat(1,[subjectnum,1]),repmat(2,[subjectnum,1])),[2,1]);
%Intensityinfo=cat(1,repmat(0,[2*subjectnum,1]),repmat(1,[2*subjectnum,1]));
%Sideinfo=repmat(cat(1,repmat(1,[subjectnum,1]),repmat(2,[subjectnum,1])),[2,1]);

[Sidecomponent,Intensitycomponent,Bothcomponent]=TCAfeatureselect(M1,Sideinfo,Intensityinfo);
% Note that it is hard to select the specific factors which encode
% both or side within a modality, thus I divided them into different
% datasets.
% [Lasercomponent,Eleccomponent]=TCAfeaturedefined(M1);
% Sidecomponent_laser=intersect(Sidecomponent,Lasercomponent);
% Sidecomponent_elec=intersect(Sidecomponent,Eleccomponent);
% Intensitycomponent_laser=intersect(Intensitycomponent,Lasercomponent);
% Intensitycomponent_elec=intersect(Intensitycomponent,Eleccomponent);
% Bothcomponent_laser=intersect(Bothcomponent,Lasercomponent);
% Bothcomponent_elec=intersect(Bothcomponent,Eleccomponent);

figure;
subplot(1,4,1)
hold on;

time_side=M1.U{2}(:,Sidecomponent);
time_side=basecorrect(mean(time_side,2),t,-0.2,0,'subtract');
plot(t,mean(time_side,2));

time_intensity=M1.U{2}(:,Intensitycomponent);
time_intensity=basecorrect(mean(time_intensity,2),t,-0.2,0,'subtract');
plot(t,mean(time_intensity,2));
time_both=M1.U{2}(:,Bothcomponent);
time_both=basecorrect(mean(time_both,2),t,-0.2,0,'subtract');
plot(t,mean(time_both,2));

title('Temporal Weights'); legend({'Side factor','Intensity factor','Both factor'});
subplot(1,4,2); topo_both=M1.U{1}(:,Bothcomponent);topoplot(mean(topo_both,2),chanlocs); caxis([0,0.2]); 
title("Both Topograph")
subplot(1,4,3);topo_intensity=M1.U{1}(:,Intensitycomponent);topoplot(mean(topo_intensity,2),chanlocs); caxis([0,0.2]);
title("Intensity Topograph");
subplot(1,4,4);
hold on;
freq=M1.U{3};
%freq=basecorrect(freq,1:6,1,6,'relativepower');
freq_Both=freq(:,Bothcomponent);
%freq_side=basecorrect(mean(freq_side,2),1:6,1,6,'normalized');
bar(1:6,mean(freq_Both,2),0.3,'r');
freq_intensity=freq(:,Intensitycomponent);
%freq_intensity=basecorrect(mean(freq_intensity,2),1:6,1,6,'normalized');
bar((1:6)+0.5,mean(freq_intensity,2),0.3,'b'); title('Freq Weights'); legend({'Both factor','Intensity factor'});


%% rank stablity for fourmodality dataset.

rankrange=[15,16,17,18,19];
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
Sideinfo=cat(1,repmat(0,[2*subjectnum,1]),repmat(1,[2*subjectnum,1]));
Intensityinfo=repmat(cat(1,repmat(1,[subjectnum,1]),repmat(2,[subjectnum,1])),[2,1]);
[Sidecomponent,Intensitycomponent,Bothcomponent]=TCAfeatureselect(M1,Sideinfo,Intensityinfo);
time_intensity=M1.U{2}(:,Intensitycomponent); 
time_both=M1.U{2}(:,Bothcomponent)
tmppeak=findpeaks(mean(time_intensity,2));
tmp=min(find(tmppeak.loc>201));
intenpeak_rank(r)=t(tmppeak.loc(tmp));
tmppeak=findpeaks(mean(time_both,2));
tmp=min(find(tmppeak.loc>201));
bothpeak_rank(r)=t(tmppeak.loc(tmp));
freq=M1.U{3};
freq_both_rank(:,r)=mean(freq(:,Bothcomponent),2);
freq_inten_rank(:,r)=mean(freq(:,Intensitycomponent),2);
end
figure;plot(rankrange, intenpeak_rank); hold on; plot(rankrange,bothpeak_rank);
figure;surf(freq_inten_rank,repmat(0.5,[size(freq_inten_rank),3]));
hold on;surf(freq_both_rank,repmat(0.8,[size(freq_inten_rank),3]));

%% summarized the factors from tensorpy_resample_laser.h5 and plot the permutated temporal factors which encode intensity or side
t=linspace(-0.2,0.5,701);
neuronsidecomptmp=[];tempsidecomptmp=[];
neuronintencomptmp=[];tempintencomptmp=[];
Sideinfo=cat(1,repmat(1,[150,1]),repmat(0,[150,1]));
Intensityinfo=repmat(cat(1,repmat(1,[75,1]),repmat(2,[75,1])),[2,1]);
for i=1:20 % 
    M1=tensor_py_resample{i,1};
    M1=normalize(M1);
    % for laser modality, there exist a light induced spike activities, remove
    % the max peak occured before 0.1s
    [Sidecomponent,Intensitycomponent,Bothcomponent]=TCAfeatureselect(M1,Sideinfo,Intensityinfo); 
    time_side=M1.U{2}(:,Sidecomponent);
     time_inten=M1.U{2}(:,Intensitycomponent);
     time_both=M1.U{2}(:,Bothcomponent);
    %
 
     timeintencomptmp(:,i)=mean(time_inten,2);
     %timeintencomptmp(:,i)=basecorrect(timeintencomptmp(:,i),t,-0.2,0,'subtract');
    %time_side=basecorrect(time_side,t,-0.2,0.1,'subtract');
    %timesidecomptmp(:,i)=mean(time_side,2);
    
    %time_both=basecorrect(mean(time_both,2),t,-0.2,0,'subtract');
    timebothcomptmp(:,i)=mean(time_both,2);
    
    freq_both=M1.U{3}(:,Bothcomponent);
    freq_side=M1.U{3}(:,Sidecomponent);
    freq_inten=M1.U{3}(:,Intensitycomponent);
    freqintencomptmp(:,i)=mean(freq_inten,2);
    freqsidecomptmp(:,i)=mean(freq_side,2);
    freqbothcomptmp(:,i)=mean(freq_both,2);
end
timeintencomptmp=basecorrect(timeintencomptmp,t,-0.2,0,'subtract');
timebothcomptmp=basecorrect(timebothcomptmp,t,-0.2,0,'subtract');
figure; subplot(1,2,1);
plot(t',smooth(nanmean(timeintencomptmp,2)),'r');hold on;
plot(t',smooth(nanmean(timebothcomptmp,2)),'g');hold on;
% plot(t',smooth(nanmean(timesidecomptmp,2)),'b');hold on;

shadebar(t',timeintencomptmp,'r');hold on;
shadebar(t',timebothcomptmp,'g');
xlim([-0.2,0.5]);
subplot(1,2,2)
plot(1:6,nanmean(freqintencomptmp,2),'r');hold on;
plot(1:6,nanmean(freqbothcomptmp,2),'g');hold on;
shadebar(1:6,freqintencomptmp,'r');hold on;
shadebar(1:6,freqbothcomptmp,'g');

%%
% t=tensormat.spkt;
tpeak=find(t>0);
invalid=isnan(timeintencomptmp(1,:))|isnan(timebothcomptmp(1,:));
timeintencomptmp(:,invalid)=[];
timebothcomptmp(:,invalid)=[];
for i=1:size(timebothcomptmp,2)
    tmppeak=findpeaks(timeintencomptmp(:,i));
    tmp=min(find(tmppeak.loc>min(tpeak)));
    intenpeak(:,i)=t(tmppeak.loc(tmp));
    tmppeak=findpeaks(timebothcomptmp(:,i));
    tmp=min(find(tmppeak.loc>min(tpeak)));
    bothpeak(:,i)=t(tmppeak.loc(tmp));
end
[~,p]=ttest(intenpeak,bothpeak);


function [Lasercomponent,Eleccomponent]=TCAfeaturedefined(M1)
 % using the accurarcy to classifiy the modalities of TCs
%     trialcomponent=M1.U{3};
    M1=normalize(M1,[],2);
     Lasercomponent=[];Eleccomponent=[];
%     truelabel=cat(1,ones(60,1),2*ones(60,1));
%     for i=1:size(trialcomponent,2)
%         mdl=fitcsvm(trialcomponent(:,i),cat(1,ones(60,1),2*ones(60,1)),'KFold',5);
%         predictlabel=kfoldPredict(mdl);
%         v=confusionmat(truelabel,predictlabel);
%         v=diag(v);
%         if v(1)>60*0.85&&v(2)<60*0.85
%             Lasercomponent=cat(1,Lasercomponent,i);
%         end
%         if v(2)>60*0.85&&v(1)<60*0.85
%             Eleccomponent=cat(1,Eleccomponent,i);
%         end
%     end
    
    modalitycomponent=M1.U{5};
    modalitycomponent=basecorrect(modalitycomponent,1:2,1,2,'relativepower');
    for i=1:size(modalitycomponent,2)
        if modalitycomponent(1,i)>0.85
            Lasercomponent=cat(1,Lasercomponent,i);
        end
        if modalitycomponent(2,i)>0.85
            Eleccomponent=cat(1,Eleccomponent,i);
        end
%          if modalitycomponent(3,i)>0.85
%             Visualcomponent=cat(1,Visualcomponent,i);
%          end
%          if modalitycomponent(4,i)>0.85
%             Auditorycomponent=cat(1,Auditorycomponent,i);
%         end
    end
            
end
function [Side,Intensity,Both]=TCAfeatureselect(M1,Sideinfo,Intensityinfo)
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
statp=reshape(statp,size(M1.U{3},2),[]);
indexvalid=[];
for i=1:2
%     if i==1
%         indexvalid(:,i)=statp(:,1)<0.05&statp(:,2)>0.05;
%     elseif i==2
%         indexvalid(:,i)=statp(:,2)<0.05&statp(:,1)>0.05;
%     elseif i==3
%         indexvalid(:,i)=statp(:,1)<0.05&statp(:,2)<0.05;
%     end
   indexvalid(:,i)=statp(:,i)<0.01;
end
Side=find(indexvalid(:,1)==1&indexvalid(:,2)~=1);
Intensity=find(indexvalid(:,2)==1&indexvalid(:,1)~=1);
Both=find(indexvalid(:,1)==1&indexvalid(:,2)==1);
end

            


