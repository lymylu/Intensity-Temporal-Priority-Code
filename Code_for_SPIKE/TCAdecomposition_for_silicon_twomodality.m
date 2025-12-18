%%
% THE tensormatrix_initial is the matrixs of channel*time*freq*(trial*subject), rank is the number of components
%tensormatrix_initial=(tensormatrix_initial-repmat(min(reshape(tensormatrix_initial,[],1,size(tensormatrix_initial,3)),[],1),[size(tensormatrix_initial,1),size(tensormatrix_initial,2),1]))./...
%    (repmat((max(reshape(tensormatrix_initial,[],1,size(tensormatrix_initial,3)),[],1)-min(reshape(tensormatrix_initial,[],1,size(tensormatrix_initial,3)),[],1)),[size(tensormatrix_initial,1),size(tensormatrix_initial,2),1]));
% % tensormatrix were normalized across each neuron.
%% transfer the tensor.h5 result from TCA_cal.py to ktensor
tensormat=matfile('/mnt/Share/yuelp/TCAdecomposition/SPKdata.mat');
tensormatrix_cal=tensormat.spike_normalized;
tensormatrix_cal=tensormatrix_cal(:,:,11:130);
trialtypelabel=tensormat.trialtypelabel;
trialtypelabel=trialtypelabel(11:130);
tensormatrix_cal1(:,:,:,1)=tensormatrix_cal(:,:,1:60);
tensormatrix_cal1(:,:,:,2)=tensormatrix_cal(:,:,61:120);
tensormatrix_cal=tensor(tensormatrix_cal1);
tensorfilename=['/mnt/Share/yuelp/TCAdecomposition/tensor_silicon_laser_elec.h5'];
tensor_py=[];rmse_py=[];similarity_py=[];
for i=1:30
    lambda=h5read(tensorfilename,['/factor_lam_',num2str(i)]);
    factortime=h5read(tensorfilename,['/factor_time_',num2str(i)]);
    factortrial=h5read(tensorfilename,['/factor_trial_',num2str(i)]);
    factorneuron=h5read(tensorfilename,['/factor_neuron_',num2str(i)]);
    factormodality=h5read(tensorfilename,['/factor_modality_',num2str(i)]);
    rmse_tmp=h5read(tensorfilename,['/reconerror_',num2str(i)]);
    rmse_py(i,:)=rmse_tmp;
    similarity_tmp=h5read(tensorfilename,['/similarity_',num2str(i)]);
    similarity_py(i,:)=similarity_tmp;
    for c=1:10
        tensor_py{i,c}=ktensor(squeeze(lambda(:,c)),squeeze(factorneuron(:,:,c))',squeeze(factortime(:,:,c))',squeeze(factortrial(:,:,c))',squeeze(factormodality(:,:,c))');
    end
end
figure; subplot(2,2,1);scatter(1:30,rmse_py); hold on; plot(1:30,mean(rmse_py,2)); title('reconstructed error of different TCA model rank');
subplot(2,2,2);scatter(2:30,diff(rmse_py)); hold on; plot(2:30,mean(diff(rmse_py),2)); title('Difference of the reconstructed error');
subplot(2,2,3);scatter(1:30,similarity_py); hold on; plot(1:30,mean(similarity_py(:,2:10),2)); title('similarity of different TCA model rank');
% svm decode accuarray using Tensor_all to identify the best rank
eventlistall=trialtypelabel(1:60);
for R=1:30
   for iter=1:10
       mdl=fitcecoc(tensor_py{R,iter}.U{3},eventlistall,'Coding','onevsall');
       label=predict(mdl,tensor_py{R,iter}.U{3}); 
       accuarray(R,iter)=sum(cellfun(@(x,y) strcmp(x,y),label,eventlistall,'UniformOutput',1))/120;
    end
end
subplot(2,2,4); hold off;scatter(1:30,squeeze(mean(accuarray,2))); title('SVM predict performance of all 4 conditions');
for i=1:10
v=locfit((1:30)',accuarray(:,i));
accuarraynew(:,i)=predict(v,(1:30)');
end
hold on; plot(1:30,mean(accuarraynew,2));
%idx_of_result=knee_pt(accuarraynew(1:end,1)); 
%%
%  generate the sampled dataspike matrix (trials were resampled)
tensormat=matfile('/mnt/Share/yuelp/TCAdecomposition/SPKdata.mat','Writable',true);
tensormatrix=tensormat.spike_normalized_laser_elec;
trialnum=[10,15,15,15,15,15,15,15,15];
data_resample=[];index_all=[];
for i=1:20
    group1index=randperm(15,10);
    data_resample(:,:,1:10,:,i)=tensormatrix(:,:,group1index,:);
    group1index=randperm(15,10);
    data_resample(:,:,11:20,:,i)=tensormatrix(:,:,group1index+15,:);
    group1index=randperm(15,10);
    data_resample(:,:,21:30,:,i)=tensormatrix(:,:,group1index+30,:);
    group1index=randperm(15,10);
    data_resample(:,:,31:40,:,i)=tensormatrix(:,:,group1index+45,:);
end
tensormat.spike_normalized_resample_laser_elec=data_resample;
%% transfer the tensorpy_resample.h5 to ktensor
tensorfilename=['/mnt/Share/yuelp/TCAdecomposition/tensor_silicon_resample_laser_elec.h5'];
tensor_py_resample=[];
for i=1:20
    lambda=h5read(tensorfilename,['/factor_lam_',num2str(i-1)]);
    factortime=h5read(tensorfilename,['/factor_time_',num2str(i-1)]);
    factortrial=h5read(tensorfilename,['/factor_trial_',num2str(i-1)]);
    factorneuron=h5read(tensorfilename,['/factor_neuron_',num2str(i-1)]);
    factormodality=h5read(tensorfilename,['/factor_modality_',num2str(i-1)]);
    for c=1:10
        tensor_py_resample{i,c}=ktensor(squeeze(lambda(:,c)),squeeze(factorneuron(:,:,c))',squeeze(factortime(:,:,c))',squeeze(factortrial(:,:,c))',squeeze(factormodality(:,:,c))');
    end
 end
%%
vizopts = {'PlotCommands',{'bar','line','scatter','bar'},...
'ModeTitles',{'Neuron','Time','Trialtype','modality'},...
'BottomSpace',0.10,'HorzSpace',0.04,'YLim',{'same','same','same','same'},'YTicks',false,'HorzSpace',0.02};
info1 = viz(normalize(tensor_py{11,1}),'Figure',1,vizopts{:});
SPKinformation=tensormat.SPKinfomation;
neurontypelabel=SPKinformation.putativeCellType;
layertypelabel=SPKinformation.Layer;
for i=1:length(info1.FactorAxes)
     if i==length(info1.FactorAxes)
         axes(info1.FactorAxes(1,i));
         replot_NeuronType(gca,neurontypelabel,layertypelabel,1);
         axes(info1.FactorAxes(3,i));
         replot_TrialType(gca,trialtypelabel,1);
     else
         axes(info1.FactorAxes(1,i));
         replot_NeuronType(gca,neurontypelabel,layertypelabel,0);
         axes(info1.FactorAxes(3,i));
         replot_TrialType(gca,trialtypelabel,0);
     end
end
%% plot the Laser induced Time TCs. Rank=11
% 1. select the modality specific TCs 
% 2. plot the temperal factor which encode side and intensity.
t=tensormat.spkt;
t=t(t>=-0.1&t<=0.3);
M1=tensor_py{11,1};
%M1=normalize(M1,[],2);%t=linspace(-0.1,0.3,400);

M1=normalize(M1);
%info1 = viz(M1,'Figure',1,vizopts{:});
[Lasercomponent,Eleccomponent]=TCAfeaturedefined(M1);
% Mlaser=extract(M1,Lasercomponent); 
timepy=M1.U{2}; validindex=[];peakindex=[];
for i=1:size(timepy,2)
    [~,index]=max(timepy(:,i));
    peakindex(i)=t(index);
end
% for laser modality, there may exist a light(laser) induced spike activities, remove
% the max peak occured before 0.1s
invalidindex=find(peakindex<0.1);
Lasercomponent(ismember(Lasercomponent,invalidindex))=[];
% for electric modality, the pulse induced artifact should be removed (the max peak occurred before 0.003s).
invalidindex=find(peakindex<0.003);
Eleccomponent(ismember(Eleccomponent,invalidindex))=[];

Intensityinfo=cat(1,repmat({'High'},[30,1]),repmat({'Low'},[30,1]));
Intensityinfo=cat(1,repmat(1,[30,1]),repmat(0,[30,1]));
Sideinfo=repmat(cat(1,repmat({'Left'},[15,1]),repmat({'Right'},[15,1])),[2,1]);
Sideinfo=repmat(cat(1,repmat(1,[15,1]),repmat(0,[15,1])),[2,1]);
[Sidecomponent,Intensitycomponent,Bothcomponent]=TCAfeatureselect(M1,Sideinfo,Intensityinfo);
Sidecomponent_laser=intersect(Sidecomponent,Lasercomponent);
Sidecomponent_elec=intersect(Sidecomponent,Eleccomponent);
Intensitycomponent_laser=intersect(Intensitycomponent,Lasercomponent);
Intensitycomponent_elec=intersect(Intensitycomponent,Eleccomponent);
Bothcomponent_laser=intersect(Bothcomponent,Lasercomponent);
Bothcomponent_elec=intersect(Bothcomponent,Eleccomponent);
figure; 
subplot(1,2,1)
hold on;
time_intensity=M1.U{2}(:,Intensitycomponent_laser);
time_intensity=basecorrect(time_intensity,t,-0.1,-0.02,'subtract');
try
    plot(t,mean(time_intensity,2));
end
time_both=M1.U{2}(:,Bothcomponent_laser);
time_both=basecorrect(time_both,t,-0.1,-0.0,'subtract');
plot(t,mean(time_both,2));
time_side=M1.U{2}(:,Sidecomponent_laser);
time_side=basecorrect(time_side,t,-0.1,-0.0,'subtract');
try plot(t,mean(time_side,2)); end
title('Temporal Weights'); legend({'Intensity factor','Both factor','Side factor'});
% time_side=double(full(extract(M1,Sidecomponent)));plot(t,squeeze(mean(mean(time_side,1),3)));
% time_intensity=double(full(extract(M1,Intensitycomponent)));plot(t,squeeze(mean(mean(time_intensity,1),3)));
subplot(1,2,2)
hold on;

time_intensity=M1.U{2}(:,Intensitycomponent_elec);
time_intensity=basecorrect(mean(time_intensity,2),t,-0.1,-0.02,'subtract');
plot(t,mean(time_intensity,2));
time_both=M1.U{2}(:,Bothcomponent_elec);
time_both=basecorrect(mean(time_both,2),t,-0.1,-0.02,'subtract');
plot(t,mean(time_both,2));
time_side=M1.U{2}(:,Sidecomponent_elec);
time_side=basecorrect(mean(time_side,2),t,-0.1,-0.02,'subtract');
plot(t,mean(time_side,2));
title('Temporal Weights'); legend({'Intensity factor','both factor','Side factor'});
%% description of neuron factors of laser and electrical stimulation.
figure; 
subplot(1,2,1);
neuron_both_laser=M1.U{1}(:,Bothcomponent_laser);
neuron_inten_laser=M1.U{1}(:,Intensitycomponent_laser);
bar(mean(neuron_both_laser,2)); title('Both laser');
replot_NeuronType(gca,neurontypelabel,layertypelabel,1);
subplot(1,2,2);
bar(mean(neuron_inten_laser,2));  title('Intensity laser');
replot_NeuronType(gca,neurontypelabel,layertypelabel,1);
%index=neuron_both_laser>prctile(neuron_both_laser,0)&neuron_inten_laser>prctile(neuron_inten_laser,0);
%[p,tbl]=anovan(neuron_inten_laser(index),{neurontypelabel(index);layertypelabel(index)},'model','full');
% [p,tbl]=anovan(neuron_both_laser,{neurontypelabel;layertypelabel},'model','full');
neurontypelabel2=double(categorical(neurontypelabel));
layertypelabel2=double(categorical(layertypelabel));
invalid=layertypelabel2==1;
[table_laser,rm_laser]=simple_mixed_anova(cat(2,neuron_both_laser(~invalid),neuron_inten_laser(~invalid)),[neurontypelabel2(~invalid);layertypelabel2(~invalid)]',{'factor'},{'neurontype','layer'});
neurotypebyfactor_laser=multcompare(rm_laser,'neurontype','By','factor','ComparisonType','bonferroni');
layerbyfactor_laser=multcompare(rm_laser,'layer','By','factor','ComparisonType','bonferroni');
factor_laserbyneuron=multcompare(rm_laser,'factor','By','neurontype','ComparisonType','bonferroni');
factor_laserbylayer=multcompare(rm_laser,'factor','By','layer','ComparisonType','bonferroni');
figure; 
neuron_both_elec=mean(M1.U{1}(:,Bothcomponent_elec),2); 
neuron_side_elec=mean(M1.U{1}(:,Sidecomponent_elec),2);

[table_elec,rm_elec]=simple_mixed_anova(cat(2,neuron_both_elec(~invalid),neuron_side_elec(~invalid)),[neurontypelabel2(~invalid);layertypelabel2(~invalid)]',{'factor'},{'neurontype','layer'});
layerbyfactor_elec=multcompare(rm_elec,'layer','By','factor','ComparisonType','bonferroni');
factor_elecbylayer=multcompare(rm_elec,'factor','By','layer','ComparisonType','bonferroni');
subplot(1,2,1);
bar(mean(neuron_both_elec,2));title('Both elec');
replot_NeuronType(gca,neurontypelabel,layertypelabel,1);
subplot(1,2,2);
bar(mean(neuron_side_elec,2));title('Side elec');
replot_NeuronType(gca,neurontypelabel,layertypelabel,1);
layertype=tabulate(layertypelabel);
layertype=layertype(:,1);
for i=2:length(layertype)
    index=~invalid&layertypelabel2==i;
    [~,p_layer(i)]=ttest(neuron_both_laser(index),neuron_inten_laser(index));
end
for i=1:length(layertype)
    index=ismember(layertypelabel,layertype{i});
    [r_laser(1,i),p_laser(1,i)]=corr(neuron_both_laser(index),neuron_inten_laser(index));
    [r_elec(1,i),p_elec(1,i)]=corr(neuron_both_elec(index),neuron_side_elec(index));
end
p_laser=mafdr(p_laser,'BHFDR',true);
p_elec=reshape(mafdr(p_elec(:),'BHFDR',true),size(p_elec));
figure; subplot(2,2,1); heatmap(r_laser);subplot(2,2,2); heatmap(p_laser);
subplot(2,2,3); heatmap(r_elec); subplot(2,2,4); heatmap(p_elec);
%%
neuronweights_laser_inten=averaged_neuronweights(neuron_inten_laser,neurontypelabel,layertypelabel);
neuronweights_laser_both=averaged_neuronweights(neuron_both_laser,neurontypelabel,layertypelabel);
%neuronweights_elec_inten=averaged_neuronweights(neuron_inten_elec,neurontypelabel,layertypelabel);
neuronweights_elec_both=averaged_neuronweights(neuron_both_elec,neurontypelabel,layertypelabel);
neuronweights_elec_side=averaged_neuronweights(neuron_side_elec,neurontypelabel,layertypelabel);
%% feature stability across rank
rankrange=[9,10,11,12,13];
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


%% pearson correlations across differenct layers and neuron types
Layer=unique(layertypelabel,'stable');
neurontype=unique(neurontypelabel,'stable');
for j=1:length(neurontype)
    for i=1:length(Layer)
        index=ismember(layertypelabel,Layer{i})&ismember(neurontypelabel,neurontype{j});
        [r(i,j),p(i,j)]=corr(neuron_inten_laser(index),neuron_inten_elec(index));
    end
end
figure; subplot(1,2,1);heatmap(r);subplot(1,2,2); heatmap(p);


%%
figure;
for i=1:2
    subplot(1,3,i);
    bar(mean(neuroncomponent{i},2));
    replot_NeuronType(gca,neurontypelist);Laser
end

%% summarized the factors from tensorpy_resample_laser.h5 and plot the permutated temporal factors which encode intensity or side
t=linspace(-0.1,0.3,400);    
neuronsidecomptmp=[];tempsidecomptmp=[];
neuronintencomptmp=[];tempintencomptmp=[];
Intensityinfo=cat(1,repmat({'High'},[20,1]),repmat({'Low'},[20,1]));
Intensityinfo=cat(1,repmat(1,[20,1]),repmat(0,[20,1]));
Sideinfo=repmat(cat(1,repmat({'Left'},[10,1]),repmat({'Right'},[10,1])),[2,1]);
Sideinfo=repmat(cat(1,repmat(1,[10,1]),repmat(0,[10,1])),[2,1]);
for i=1:20 % 
    M1=tensor_py_resample{i,1};
    M1=normalize(M1);
    % for laser modality, there exist a light induced spike activities, remove
    % the max peak occured before 0.1s
    [Lasercomponent,Eleccomponent]=TCAfeaturedefined(M1);
    % Mlaser=extract(M1,Lasercomponent); 
    timepy=M1.U{2}; validindex=[];peakindex=[];
    for c=1:size(timepy,2)
        [~,index]=max(timepy(:,c));
        peakindex(c)=t(index);
    end
    invalidindex=find(peakindex<0.1);
    Lasercomponent(ismember(Lasercomponent,invalidindex))=[];
    invalidindex=find(peakindex<0.005);
    Eleccomponent(ismember(Eleccomponent,invalidindex))=[];
    [Sidecomponent,Intensitycomponent,Bothcomponent]=TCAfeatureselect(M1,Sideinfo,Intensityinfo); 
    time_side_laser=M1.U{2}(:,intersect(Sidecomponent,Lasercomponent));
    time_side_elec=M1.U{2}(:,intersect(Sidecomponent,Eleccomponent));
    time_side_laser=basecorrect(mean(time_side_laser,2),t,-0.1,-0.02,'subtract');
    time_side_elec=basecorrect(mean(time_side_elec,2),t,-0.1,-0.02,'subtract');
     time_inten_laser=M1.U{2}(:,intersect(Intensitycomponent,Lasercomponent));
    time_inten_elec=M1.U{2}(:,intersect(Intensitycomponent,Eleccomponent));
    time_inten_laser=basecorrect(mean(time_inten_laser,2),t,-0.1,-0.02,'subtract');
    time_inten_elec=basecorrect(mean(time_inten_elec,2),t,-0.1,-0.02,'subtract');
     time_both_laser=M1.U{2}(:,intersect(Bothcomponent,Lasercomponent));
    time_both_elec=M1.U{2}(:,intersect(Bothcomponent,Eleccomponent));
    time_both_laser=basecorrect(mean(time_both_laser,2),t,-0,-0,'subtract');
    time_both_elec=basecorrect(mean(time_both_elec,2),t,-0,-0,'subtract');
    timeintencomptmp_laser(:,i)=mean(time_inten_laser,2);
   try
    timesidecomptmp_laser(:,i)=mean(time_side_laser,2);
   end
    timebothcomptmp_laser(:,i)=mean(time_both_laser,2);
    try
    timeintencomptmp_elec(:,i)=mean(time_inten_elec,2);
    end
    timesidecomptmp_elec(:,i)=mean(time_side_elec,2);
    timebothcomptmp_elec(:,i)=mean(time_both_elec,2);
    neuronsidecomp_laser=M1.U{1}(:,intersect(Sidecomponent,Lasercomponent));
    neuronsidecomp_elec=M1.U{1}(:,intersect(Sidecomponent,Eleccomponent));
    neuronintencomp_laser=M1.U{1}(:,intersect(Intensitycomponent,Lasercomponent));
    neuronintencomp_elec=M1.U{1}(:,intersect(Intensitycomponent,Eleccomponent));
    neuronbothcomp_laser=M1.U{1}(:,intersect(Bothcomponent,Lasercomponent));
    neuronbothcomp_elec=M1.U{1}(:,intersect(Bothcomponent,Eleccomponent));
    neurontypelabel=repmat({'all'},[1,235]);
    %neuronweights_laser_inten(i,:,:)=averaged_neuronweights(neuronintencomp_laser,neurontypelabel,layertypelabel);
    %neuronweights_laser_both(i,:,:)=averaged_neuronweights(neuronbothcomp_laser,neurontypelabel,layertypelabel);
    %neuronweights_elec_inten(i,:,:)=averaged_neuronweights(neuronintencomp_elec,neurontypelabel,layertypelabel);
    %neuronweights_elec_both(i,:,:)=averaged_neuronweights(neuronbothcomp_elec,neurontypelabel,layertypelabel);
    %neuronweights_elec_side(i,:,:)=averaged_neuronweights(neuronsidecomp_elec,neurontypelabel,layertypelabel);
end
figure; subplot(1,2,1);
plot(t',smooth(nanmean(timeintencomptmp_laser,2)),'r');hold on;
plot(t',smooth(nanmean(timebothcomptmp_laser,2)),'g');hold on;
shadebar(t',timeintencomptmp_laser,'r');hold on;
shadebar(t',timebothcomptmp_laser,'g');
subplot(1,2,2);
plot(t',smooth(nanmean(timeintencomptmp_elec,2)),'r');hold on;
plot(t',smooth(nanmean(timesidecomptmp_elec,2)),'b');hold on;
plot(t',smooth(nanmean(timebothcomptmp_elec,2)),'g');hold on;


shadebar(t',timeintencomptmp_elec,'r');hold on;
shadebar(t',timesidecomptmp_elec,'b');
shadebar(t',timebothcomptmp_elec,'g');
xlim([-0.05,0.1]);

[table_elec,rm_elec]=simple_mixed_anova(cat(2,neuronweights_elec_both,neuronweights_elec_side),[],{'factor','layer'});
layerbyfactor_elec=multcompare(rm_elec,'layer','By','factor','ComparisonType','bonferroni');
factor_elecbylayer=multcompare(rm_elec,'factor','By','layer','ComparisonType','bonferroni');

[table_laser,rm_laser]=simple_mixed_anova(cat(2,neuronweights_laser_both,neuronweights_laser_inten),[],{'factor','layer'});
layerbyfactor_laser=multcompare(rm_laser,'layer','By','factor','ComparisonType','bonferroni');
factor_laserbylayer=multcompare(rm_laser,'factor','By','layer','ComparisonType','bonferroni');



%%
% t=tensormat.spkt;
tpeak=find(t>0.11);

invalid=isnan(timeintencomptmp_laser(1,:))|isnan(timebothcomptmp_laser(1,:));
timeintencomptmp_laser(:,invalid)=[];
timebothcomptmp_laser(:,invalid)=[];
for i=1:size(timebothcomptmp_laser,2)
    tmppeak=findpeaks(timeintencomptmp_laser(:,i));
    tmp=min(find(tmppeak.loc>=min(tpeak)));
    intenpeak_laser(:,i)=t(tmppeak.loc(tmp));
    tmppeak=findpeaks(timebothcomptmp_laser(:,i));
    tmp=min(find(tmppeak.loc>=min(tpeak)));
    bothpeak_laser(:,i)=t(tmppeak.loc(tmp));
end
[~,p_laser]=ttest(intenpeak_laser,bothpeak_laser);
% figure;
% for i=1:size(timebothcomptmp_laser,2)
%     subplot(4,5,i)
%     plot(t,timeintencomptmp_laser(:,i),'r');hold on; plot(t,timebothcomptmp_laser(:,i),'b');
% end
tpeak=find(t>=0);
invalid=isnan(timeintencomptmp_elec(1,:))|isnan(timebothcomptmp_elec(1,:))|isnan(timesidecomptmp_elec(1,:));
timeintencomptmp_elec(:,invalid)=[];
timebothcomptmp_elec(:,invalid)=[];
timesidecomptmp_elec(:,invalid)=[];
for i=1:size(timebothcomptmp_elec,2)
    tmppeak=findpeaks(timeintencomptmp_elec(:,i));
    tmp=min(find(tmppeak.loc>=min(tpeak)));
    intenpeak_elec(:,i)=t(tmppeak.loc(tmp));
    tmppeak=findpeaks(timebothcomptmp_elec(:,i));
    tmp=min(find(tmppeak.loc>=min(tpeak)));
    bothpeak_elec(:,i)=t(tmppeak.loc(tmp));
    tmppeak=findpeaks(timesidecomptmp_elec(:,i));
    tmp=min(find(tmppeak.loc>=min(tpeak)));
    sidepeak_elec(:,i)=t(tmppeak.loc(tmp));
end
[~,p_elec]=ttest(intenpeak_elec,bothpeak_elec);


function [Lasercomponent,Eleccomponent]=TCAfeaturedefined(M1)
 % using the accurarcy to classifiy the modalities of TCs
%     trialcomponent=M1.U{3};
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
    
    modalitycomponent=M1.U{4};
    modalitycomponent=basecorrect(modalitycomponent,1:2,1,2,'relativepower');
    for i=1:size(modalitycomponent,2)
        if modalitycomponent(1,i)>0.85
            Lasercomponent=cat(1,Lasercomponent,i);
        end
        if modalitycomponent(2,i)>0.85
            Eleccomponent=cat(1,Eleccomponent,i);
        end
    end
            
end
function [Side,Intensity,Both]=TCAfeatureselect(M1,Sideinfo,Intensityinfo)
% return the component index which significantly related to Side or Intensity.
    trialcomponent=M1.U{3}; % using trial factor to fit.
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
function replot_TrialType(figaxes,trialtypelabel,xtickopt)
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
 a=tabulate(trialtypelabel);
 bounddepth=[0;cumsum(cell2mat(a(:,2)))];
 hold on;
 for j=1:length(bounddepth)-1
    if j>1
    plot([bounddepth(j),bounddepth(j)],figaxes.YLim,'black','LineStyle','--');
    end
 end
 if xtickopt==1
 xtick=arrayfun(@(x,y) (y-x)/2+x,bounddepth(1:end-1),bounddepth(2:end),'UniformOutput',1);
 set(gca,'XTick',xtick,'XTickLabel',trialtypes);
 end
end
function replot_NeuronType(figaxes,celltypelabel,layertypelabel,xtickopt)
% replot the current figure to differentiate cell types and layer positions from the neuron
% associate the bar plot of TCs and reconstructed TCs
 axes(figaxes);
 data=findobj(figaxes);
 Y_data=data(end).YData;
 X_data=data(end).XData;
 delete(data(end));hold on;
 bar(X_data(ismember(celltypelabel,'Pyramidal Cell')),Y_data(ismember(celltypelabel,'Pyramidal Cell')),'b','BarWidth',0.5,'EdgeAlpha',0);
 bar(X_data(ismember(celltypelabel,'Interneuron')),Y_data(ismember(celltypelabel,'Interneuron')),'r','BarWidth',0.5,'EdgeAlpha',0);
 a=tabulate(layertypelabel);
 bounddepth=[0;cumsum(cell2mat(a(:,2)))];
 hold on;
 for j=1:length(bounddepth)-1
    if j>1
    plot([bounddepth(j),bounddepth(j)],figaxes.YLim,'black','LineStyle','--');
    end
 end
 if xtickopt==1
 xtick=arrayfun(@(x,y) (y-x)/2+x,bounddepth(1:end-1),bounddepth(2:end),'UniformOutput',1);
 set(gca,'XTick',xtick,'XTickLabel',{'I','II/III','IV','Va','Vb','VI'});
 end
end
function neuronweights=averaged_neuronweights(TC,neurontype,layertype)

neuronsubtype=unique(neurontype);
layersubtype=unique(layertype);
for i=1:length(neuronsubtype)
    for j=1:length(layersubtype)
        try
        neuronweights{i,j}=TC(ismember(neurontype,neuronsubtype{i})&ismember(layertype,layersubtype{j}));
        catch
            neuronweights(i,j)=nan;
        end
        
        end
end  
end


