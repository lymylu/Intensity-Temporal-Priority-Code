%% main function of the S1 silicon SPK data
cd '/media/huli/MyBook/YLP/siliconS1/'
datamat=matfile('/media/huli/MyBook/YLP/Neuraltrajectory/workspaceofReAnalysis.mat','Writable',true);
subjectnum={'rat7silicon','rat9silicon','rat10silicon','rat10silicon2','rat11silicon','rat12silicon','rat13silicon','rat13silicon2','rat14silicon','rat15silicon','rat15silicon2','rat15silicon3','rat16silicon','rat16silicon2','rat17silicon'}; % electrode 1-3 in rat9 are broken. padding with NAN in CSD average.
%subjectnum={'rat17silicon'};
depth=linspace(100,1600,16);
condition_laser={'Left_3_5_S1_3.mat','Right_3_5_S1_3.mat','Left_3_S1_3.mat','Right_3_S1_3.mat'};
condition_elec={'Left_1_S1_3.mat','Right_1_S1_3.mat','Left_75_S1_3.mat','Right_75_S1_3.mat'};
%spkt=linspace(-2,4,6/0.01+1); % for the bin size of 10 ms
%% loading the laser and sound data from /SPKlasernew rawdata(can be change to bindata)
% rawtopsth-noaverage,0.05 for neural trajectory, rawtobin-noaverage,0.01 for anova
for i=1:length(condition_laser)
tmp=cellfun(@(x) Loaddata('/media/huli/MyBook/YLP/siliconS1/SPKlaser_alltrial/',condition_laser{i},x,'rawtopsth-noaverage',0.005),subjectnum,'UniformOutput',0); % rawtobin or bindata
tmp=cellfun(@(x,y) LoadBehaviorScore('laserdescription.xlsx',condition_laser{i},x,y),subjectnum,tmp,'UniformOutput',0);
tmp=catall(tmp);
try
   spkt=tmp.spkt;
end
%  tmp.dataoutput=cellfun(@(x) basecorrect(x,spkt,-0.2,0,'Zscore'),tmp.dataoutput,'UniformOutput',0); %  basecorrect from 0.5s before stimulus
eval(['SPKsummary.',condition_laser{i}(1:end-4),'=tmp.dataoutput;']);
eval(['SPKsummary.',condition_laser{i}(1:end-4),'_behaviorscore=tmp.behaviorscore;']);
end
SPKsummary.spikename=tmp.spikename;
SPKsummary.spikesubject=tmp.spikesubject;
SPKsummary.spkt_laser=spkt;
SPKsummary=Rearrange_trialorder(SPKsummary,'Left_3_5_S1_3','Left_3_S1_3');
SPKsummary=Rearrange_trialorder(SPKsummary,'Right_3_5_S1_3','Right_3_S1_3');
% loading the elec data from SPKvonfreynew raw data
for i=1:length(condition_elec)
tmp=cellfun(@(x) Loaddata('/media/huli/MyBook/YLP/siliconS1/SPKvonfrey_alltrial/',condition_elec{i},x,'rawtopsth-noaverage',0.005),subjectnum,'UniformOutput',0); % rawtobin or bindata
tmp=cellfun(@(x,y) LoadBehaviorScore('electricdescription.xlsx',condition_elec{i},x,y),subjectnum,tmp,'UniformOutput',0);
tmp=catall(tmp);
try
   spkt=tmp.spkt;
end
%tmp.dataoutput=cellfun(@(x) basecorrect(x,spkt,-0.2,0,'Zscore'),tmp.dataoutput,'UniformOutput',0); %  basecorrect from 0.5s before stimulus
eval(['SPKsummary.',condition_elec{i}(1:end-4),'=tmp.dataoutput;']);
eval(['SPKsummary.',condition_elec{i}(1:end-4),'_behaviorscore=tmp.behaviorscore;']);
eval(['SPKsummary.',condition_elec{i}(1:end-4),'_behaviorscore=cat(1,cell(length(SPKsummary.Left_3_5_S1_3)-length(SPKsummary.',condition_elec{i}(1:end-4),'),1),SPKsummary.',condition_elec{i}(1:end-4),'_behaviorscore);']);
eval(['SPKsummary.',condition_elec{i}(1:end-4),'=cat(1,cell(length(SPKsummary.Left_3_5_S1_3)-length(SPKsummary.',condition_elec{i}(1:end-4),'),1),SPKsummary.',condition_elec{i}(1:end-4),');']);
end
%SPKsummary.sound_S1_3=cat(1,cell(length(SPKsummary.Left_1_S1_3)-length(SPKsummary.sound_S1_3),1),SPKsummary.sound_S1_3);
SPKsummary=Rearrange_trialorder(SPKsummary,'Left_1_S1_3','Left_75_S1_3');
SPKsummary=Rearrange_trialorder(SPKsummary,'Right_1_S1_3','Right_75_S1_3');
SPKsummary.spkt_elec=spkt;
SPKsummary=spikeproperties(SPKsummary);
% datamat.SPKsummary_psth=SPKsummary;
% %  correct the normalized depth of spk units
% SPKsummary=datamat.SPKsummary_psth;
[~,~,b1]=xlsread('./channeldepth.xlsx');
channel=1:16;
depththick=[1,1,4,3,2,3,3,1]; 
[x2,y2]=meshgrid(1,0:sum(depththick)-1);
for i=1:length(subjectnum)
    depth=cellfun(@(x) channellength(x),b1(2+i,2:9),'UniformOutput',1);
    intSPKy=[];
     linpoint={[-1,0],[0,1],[2,5],[6,8],[9,10],[11,13],[14,16],[17,18]};
    for l=1:length(depth)
        if depth(l)~=0
            intSPKy=cat(2,intSPKy,linspace(linpoint{l}(1),linpoint{l}(2),depth(l)));
        end
    end
    [x1,y1(:,i)]=meshgrid(1,intSPKy);
end
%
Layertitle={'Outup','LayerI','LayerII/III','LayerIV','LayerVa','LayerVb','LayerVI','Outdown'};
for i=1:length(SPKsummary.spikesubject)
    index=ismember(subjectnum,SPKsummary.spikesubject{i});
    depth_norm=y1(:,index);
    SPKsummary.norm_depth(i)=depth_norm(SPKsummary.maxWaveformCh1(i));
    for j=1:length(linpoint)
        if SPKsummary.norm_depth(i)>=linpoint{j}(1)&&SPKsummary.norm_depth(i)<=linpoint{j}(2)
            SPKsummary.Layer{i}=Layertitle{j};
        end
    end
end
% using the 0.38 as the classfication of putativeneuron and interneuron
for i=1:length(SPKsummary.troughToPeak)
    if SPKsummary.troughToPeak(i)<0.38
        SPKsummary.putativeCellType{i}='Interneuron';
    else
        SPKsummary.putativeCellType{i}='Pyramidal Cell';
    end
end
datamat.SPKsummary_psth=SPKsummary;
%% cat the spike martix for tensormatrix
for c=1:2
SPK=[]; 
if c==1; condition=condition_laser;
else  condition=condition_elec;end
for i=1:length(condition) % lefthigh, righthigh, leftlow, rightlow
    tmp=eval(['SPKsummary.',condition{i}(1:end-4),';']);
    SPK_laser_tmp=[];
    for j=1:length(tmp)
    if size(tmp{j},2)>15
        tmp{j}=tmp{j}(:,randperm(size(tmp{j},2),15));
    elseif size(tmp{j},2)<15
        tmp{j}=cat(2,tmp{j},nan(size(tmp{j},1),15-size(tmp{j},2)));
    end
    SPK_tmp(:,:,j)=tmp{j};
    end
    SPK=cat(2,SPK,SPK_tmp);
end
    SPK=permute(SPK,[3,1,2]);
    SPK(isnan(SPK))=0;
    if c==1; datamat.SPK_laser_reassign=SPK;
    else datamat.SPK_elec_reassign=SPK;
    end
end
datamat.spkt=SPKsummary.spkt_elec;
%%  neural trajectory analysis
SPKsummary=datamat.SPKsummary_psth;
cm=flipud(gray(256)); cm=cm(1:128,:);
% %%
%SPKsummary=SPKsummary_bin;
%spkt=spkt+0.1;
% pca decomposition
close all; 
clear dist* angle*
spkt=SPKsummary.spkt_elec;
%compare={[1,2],[3,4],[1,3],[2,4]}; % 1,2 intensity difference 3,4 Side difference
compare2={[1,2,3],[2,1,4],[3,4,1],[4,3,2]}; % normalized perferential comparisons. 1,2,3 means 1-2,1-3/ [（1-2）+（1-3）]
pcarange=[-0.2,0.5]; % long pca range need more pcs to explain the 85% varaince.
conditioncolor={[0,0.447,0.741],[0.85,0.325,0.098],[0.929,0.694,0.125],[0.494,0.184,0.556]};
Layertitle=unique(SPKsummary.Layer);
Layertitle{end+1}='All';
spktroi=spkt(spkt>pcarange(1)&spkt<pcarange(2));
%
for lay=7 % for all neurons
    if ~strcmp(Layertitle{lay},'All')
        %valid=SPKsummary.firingRate>prctile(SPKsummary.firingRate,10)&strcmp(SPKsummary.Layer,Layertitle{lay});
        %valid=SPKsummary.firingRate>prctile(SPKsummary.firingRate,10);
        %&strcmp(SPKsummary.putativeCellType,'Pyramidal Cell');
        valid=strcmp(SPKsummary.Layer,Layertitle{lay});
    else
        valid=SPKsummary.firingRate>prctile(SPKsummary.firingRate,10);
    end
    if lay==2
        Layertitle{lay}='LayerII_III';
    end

coeff=[];
SPKsummary_normalized=normalized(SPKsummary,spkt,pcarange,condition_elec,valid);
[~,score,latent,~,explain]=pca(SPKsummary_normalized); % PCA 
 figure; try subplot(1,4,1); plot(1:20,explain(1:20)); end 
[coeff,chopSPKsummary]=chopPCAmatrix(score,SPKsummary_normalized,sum(spkt>pcarange(1)&spkt<pcarange(2)));
%subplot(1,4,2); plot(spkt(spkt>pcarange(1)&spkt<pcarange(2)),mean(chopSPKsummary,3));
for i=1:4
   % [d,z]=procrustes(coeff(:,1:3,1),coeff(:,1:3,i),'scaling',false,'reflection',false);
%     plot3(z(:,1),z(:,2),z(:,3),'color',conditioncolor{i},'LineStyle','--');hold on;
%     plot3(z(5,1),z(5,2),z(5,3),'.','MarkerSize',20,'color',conditioncolor{i});
     subplot(1,4,2); plot(spkt(spkt>pcarange(1)&spkt<pcarange(2)),mean(chopSPKsummary(:,:,i),2),'color',conditioncolor{i}); hold on;
    subplot(1,4,3);
    plot3(coeff(:,1,i),coeff(:,2,i),coeff(:,3,i),'color',conditioncolor{i},'LineStyle','--'); hold on;
    % plot3(coeff(5,1,i),coeff(5,2,i),coeff(20,3,i),'.','MarkerSize',20,'color',conditioncolor{i});
end
pcs_elec=min(find(cumsum(explain)>85));% 5 components explain more than 85% for laser
pcs_elec_sum(lay)=pcs_elec;
[distchange_elec,degree_elec,angle_elec]=cal_distangle(coeff,pcs_elec,compare2);
distchange_elec_all{lay}=distchange_elec;
degree_elec_all{lay}=degree_elec;
%averageddistchange_elec(:,lay)=squeeze(mean(distchange_elec(spktroi>0&spktroi<0.1,:),1));
subplot(1,4,4); plot(spkt(spkt>pcarange(1)&spkt<pcarange(2)),distchange_elec); 
% for pcarange=[-1,2];
%savefig(['/media/huli/My Book/YLP/Neuraltrajectory/Figures/neurotrajectory_elec_',Layertitle{lay},'_longtime.fig']);
% for pcarange=[-0.2,0.5]; used in the manuscript
%savefig(['/media/huli/My Book/YLP/Neuraltrajectory/Figures/neurotrajectory_elec_',Layertitle{lay},'.fig']);
close gcf;
figure; 
for i=1:4
    subplot(1,4,i); [~,sortindex]=sort(mean(chopSPKsummary(spktroi>0&spktroi<0.1,:,1),1));
    imagesc(spktroi,1:length(sortindex),chopSPKsummary(:,sortindex,i)');caxis([-2,5]);title(condition_elec{i});
end
colormap jet;
%savefig(['/media/huli/My Book/YLP/Neuraltrajectory/Figures/neurons_elec_generalplot_longtime.fig']);
close gcf;
coeff=[];
spkt=SPKsummary.spkt_laser;
SPKsummary_normalized=normalized(SPKsummary,spkt,pcarange,condition_laser,valid);
SPKsummary_normalized(:,isnan(SPKsummary_normalized(1,:)))=[];
[~,score,latent,~,explain]=pca(SPKsummary_normalized);
 figure; try subplot(1,4,1); plot(1:20,explain(1:20));end % 8 components  explain more than 85% for electric
[coeff,chopSPKsummary]=chopPCAmatrix(score,SPKsummary_normalized,sum(spkt>pcarange(1)&spkt<pcarange(2)));
%coeff=coeff(:,:,5:8);chopSPKsummary=chopSPKsummary(:,:,5:8);
for i=1:4
 %[d,z]=procrustes(coeff(:,1:3,1),coeff(:,1:3,i),'scaling',false,'reflection',false);
%  plot3(z(:,1),z(:,2),z(:,3),'color',conditioncolor{i});hold on;
%  plot3(z(5,1),z(5,2),z(5,3),'.','MarkerSize',20,'color',conditioncolor{i});
subplot(1,4,2); plot(spkt(spkt>pcarange(1)&spkt<pcarange(2)),mean(chopSPKsummary(:,:,i),2),'color',conditioncolor{i}); hold on;
subplot(1,4,3);
    plot3(coeff(:,1,i),coeff(:,2,i),coeff(:,3,i),'color',conditioncolor{i}); hold on;
    %plot3(coeff(20,1,i),coeff(2,2,i),coeff(5,3,i),'.','MarkerSize',20,'color',conditioncolor{i});
end 
pcs_laser=min(find(cumsum(explain)>85));
pcs_laser_sum(lay)=pcs_laser;
[distchange_laser,degree_laser,angle_laser]=cal_distangle(coeff,pcs_laser,compare2);
%averageddistchange_laser(:,lay)=squeeze(mean(distchange_laser(spktroi>0.1&spktroi<0.2,:),1));
distchange_laser_all{lay}=distchange_laser;
degree_laser_all{lay}=degree_laser;
subplot(1,4,4); plot(spkt(spkt>pcarange(1)&spkt<pcarange(2)),distchange_laser);
%hold on; plot(spkt(spkt>pcarange(1)&spkt<pcarange(2)),degree_laser,'black');
%savefig(['/media/huli/My Book/YLP/Neuraltrajectory/Figures/neurotrajectory_laser_',Layertitle{lay},'_longtime.fig']);
close gcf;
figure; 
for i=1:4
    subplot(1,4,i); [~,sortindex]=sort(mean(chopSPKsummary(spktroi>0,:,1),1));
    imagesc(spktroi,1:length(sortindex),chopSPKsummary(:,sortindex,i)');caxis([-2,5]);title(condition_laser{i});
end
colormap jet;
%savefig(['/media/huli/My Book/YLP/Neuraltrajectory/Figures/neurons_laser_generalplot_longtime.fig']);
close gcf;
%%  Generate the shuffle neuraltrajectory to estimate the distribution and calculate the significant interval.
clear *shuffle*
rng('default'); % for replicable;
for m=1:20
    spkt_elec=SPKsummary.spkt_elec;
coeff=[];
SPKsummary_normalized=normalized_shuffle(SPKsummary,spkt_elec,pcarange,condition_elec,valid);
SPKsummary_normalized(:,isnan(SPKsummary_normalized(1,:)))=[];
[~,score,latent,~,explain]=pca(SPKsummary_normalized);
explain_elec(:,m)=explain(1:20);
pcs_elec=min(find(cumsum(explain)>85));
%pcs_elec=9;
%pcs_elec_shu(lay,m)=pcs_elec;
[coeff,chopSPKsummary]=chopPCAmatrix(score,SPKsummary_normalized,sum(spkt_elec>pcarange(1)&spkt_elec<pcarange(2)));
%pcs=min(find(cumsum(explain)>85));
[distchange_shuffle_elec(m,:,:)]=cal_distangle(coeff,pcs_elec,compare2);
averageddistchange_elec_shu(m,:,lay)=squeeze(mean(distchange_shuffle_elec(m,spktroi>0&spktroi<0.1,:),2));
%
spkt_laser=SPKsummary.spkt_laser;
coeff=[];
SPKsummary_normalized=normalized_shuffle(SPKsummary,spkt_laser,pcarange,condition_laser,valid);
SPKsummary_normalized(:,isnan(SPKsummary_normalized(1,:)))=[];
[~,score,latent,~,explain]=pca(SPKsummary_normalized);
explain_laser(:,m)=explain(1:20);
[coeff,chopSPKsummary]=chopPCAmatrix(score,SPKsummary_normalized,sum(spkt_laser>pcarange(1)&spkt_laser<pcarange(2)));
pcs_laser=min(find(cumsum(explain)>85));
%pcs_laser=7;
%pcs_laser_shu(lay,m)=pcs_laser;
[distchange_shuffle_laser(m,:,:)]=cal_distangle(coeff,pcs_laser,compare2);
averageddistchange_laser_shu(m,:,lay)=squeeze(mean(distchange_shuffle_laser(m,spktroi>0.1&spktroi<0.5,:),2));
% shuffle baseline
[distchange_shuffle_elecbaseline(m,:,:)]=pca_shuffle_baseline(SPKsummary,spkt_elec,pcarange,condition_elec,valid,compare2,8);
[distchange_shuffle_laserbaseline(m,:,:)]=pca_shuffle_baseline(SPKsummary,spkt_laser,pcarange,condition_laser,valid,compare2,5);
averageddistchange_elecbaseline_shu(m,:,lay)=squeeze(mean(distchange_shuffle_elecbaseline(m,spktroi>0&spktroi<0.1,:),2));
averageddistchange_laserbaseline_shu(m,:,lay)=squeeze(mean(distchange_shuffle_laserbaseline(m,spktroi>0.1&spktroi<0.5,:),2));
end
distchange_shu_elec_all{lay}=distchange_shuffle_elec;
distchange_shu_laser_all{lay}=distchange_shuffle_laser;
figure; subplot(1,2,1); shadebar(1:20,explain_elec,'black','std');subplot(1,2,2);shadebar(1:20,explain_laser,'red','std');
%savefig(['/media/huli/My Book/YLP/Neuraltrajectory/Figures/explain_shuffle.fig']);
% plot variablity
figure;c=1;  
for i=1:4
    subplot(4,2,c)
    spkt=spkt_elec;
    imagesc(gca,spkt(spkt>pcarange(1)&spkt<pcarange(2)),-20:20,repmat(degree_elec(:,i),[1,41])');colormap(cm); hold on; caxis([0,50]); axis xy;
    shadebar(spkt(spkt>pcarange(1)&spkt<pcarange(2)),distchange_shuffle_elec(:,:,i)',conditioncolor{i},'std');hold on; ylim([-10,10]);
    h=timevaryingttest(distchange_shuffle_elec(:,:,i),distchange_shuffle_elecbaseline(:,:,i));
   
    plot(spkt(spkt>pcarange(1)&spkt<pcarange(2)),h*10);
    shadebar(spkt(spkt>pcarange(1)&spkt<pcarange(2)),distchange_shuffle_elecbaseline(:,:,i)',conditioncolor{i},'std');hold on; ylim([-10,10]);
    
    c=c+1;
    subplot(4,2,c)
    spkt=spkt_laser;
    imagesc(gca,spkt(spkt>pcarange(1)&spkt<pcarange(2)),-20:20,repmat(degree_laser(:,i),[1,41])');colormap(cm); hold on; caxis([0,50]); axis xy;
    shadebar(spkt(spkt>pcarange(1)&spkt<pcarange(2)),distchange_shuffle_laser(:,:,i)',conditioncolor{i},'std');hold on;ylim([-10,10]);
    shadebar(spkt(spkt>pcarange(1)&spkt<pcarange(2)),distchange_shuffle_laserbaseline(:,:,i)',conditioncolor{i},'std');hold on;ylim([-10,10]);
    h=timevaryingttest(distchange_shuffle_laser(:,:,i),distchange_shuffle_laserbaseline(:,:,i));
    plot(spkt(spkt>pcarange(1)&spkt<pcarange(2)),-h*10);
    
    c=c+1;
end
%savefig(['/media/huli/My Book/YLP/Neuraltrajectory/Figures/neurotrajectory_shuffle_',Layertitle{lay},'_longtime.fig']);
close gcf
end
%%
filemat=matfile('./trajectorylayerstat.mat','Writable',true);
cm=jet(5);
figure;
for j=1:4
    subplot(2,2,j);
for i=1:5
    hold on;plot(spktroi,squeeze(mean(distchange_shu_laser_all{i+1}(:,:,j),1))','color',cm(i,:));
end
for i=1:5
    hold on; 
    shadebar(spktroi,distchange_shu_laser_all{i+1}(:,:,j)',cm(i,:),'std'); ylim([-10,10]); xlim([-0.2,0.5]);
end
legend(Layertitle(2:6));
end
%savefig(gcf,'./Figures/laser_layer_shuffle.fig');
%close gcf;
figure;
for j=1:4
    subplot(2,2,j);
for i=1:5
    %hold on;plot(spktroi,squeeze(mean(distchange_shu_elec_all{i+1}(:,:,j),1))','color',cm(i,:));
end

for i=1:5
    hold on; 
    shadebar(spktroi,distchange_shu_elec_all{i+1}(:,:,j)',cm(i,:),'std'); ylim([-10,10]); xlim([-0.2,0.5]);
end
legend(Layertitle(2:6));
end
%savefig(gcf,'./Figures/elec_layer_shuffle.fig');
%close gcf;
% summarized the distance within the given time interval and compare among layers
laser_distance=cellfun(@(x) squeeze(mean(x(:,spktroi>0.1,[1,2]),2)),distchange_shu_laser_all(2:6),'UniformOutput',0);
elec_distance=cellfun(@(x) squeeze(mean(x(:,spktroi>0&spktroi<0.1,[1,3]),2)),distchange_shu_elec_all(2:6),'UniformOutput',0);
laser_tmp=[];elec_tmp=[];
for i=1:length(laser_distance)
    laser_tmp(:,:,i)=laser_distance{i};
    elec_tmp(:,:,i)=elec_distance{i};
end
filemat.laser_distance=laser_distance;
filemat.elec_distance=elec_distance;
[tbl1,rm1]=simple_mixed_anova(laser_tmp,[],{'factor','layer'});
[tbl2,rm2]=simple_mixed_anova(elec_tmp,[],{'factor','layer'});
a=multcompare(rm1,'layer','by','factor');
b=multcompare(rm2,'layer','by','factor');
% save the shuffled neural trajectories 
filemat.laser_shuffle=distchange_shu_laser_all;
filemat.elec_shuffle=distchange_shu_elec_all;
filemat.layertitle=Layertitle(2:6);
figure;
for i=2:7
subplot(2,4,i) 
plot(squeeze(mean(averageddistchange_laserbaseline_shu(:,:,i),1)));hold on;
plot(squeeze(mean(averageddistchange_laser_shu(:,:,i),1)));title(Layertitle{i});
for c=1:4
[h_shu(c)]=ttest(averageddistchange_laserbaseline_shu(:,c,i),averageddistchange_laser_shu(:,c,i));
end
plot(h_shu); ylim([-15,15]);
end

figure;
for i=2:7
subplot(2,4,i) 
plot(squeeze(mean(averageddistchange_elecbaseline_shu(:,:,i),1)));hold on;
plot(squeeze(mean(averageddistchange_elec_shu(:,:,i),1)));title(Layertitle{i});
for c=1:4
[h_shu(c)]=ttest(averageddistchange_elecbaseline_shu(:,c,i),averageddistchange_elec_shu(:,c,i));
end
plot(h_shu); ylim([-15,15]);
end

%% behavior score
[~,index]=unique(SPKsummary.spikesubject);
for i=1:4
tmp=eval(['SPKsummary.',condition_laser{i+1}(1:end-4),'_behaviorscore']);
tmp=tmp(index);
laserscore(:,i)=cellfun(@(x) mean(x),tmp,'UniformOutput',1);
end
for i=1:4
tmp=eval(['SPKsummary.',condition_elec{i}(1:end-4),'_behaviorscore']);
tmp=tmp(index);
elecscore(:,i)=cellfun(@(x) mean(x),tmp,'UniformOutput',1);
end
%% plot binned spike firing % for generate figures to TCA article
close all;
datastat=matfile('/media/huli/My Book/YLP/TCAdecomposition/Figures/FiringrateStat2.mat','Writable',true);
Layers={'All','LayerII/III','LayerIV','LayerVa','LayerVb','LayerVI'};
Savelayer={'All','LayerII_III','LayerIV','LayerVa','LayerVb','LayerVI'};
SPKsummary=datamat.SPKsummary_psth;
%SPKsummary=cal_responses(SPKsummary,spkt_laser,condition_laser,[0,0.3],'Laserresponse');
%SPKsummary=cal_responses(SPKsummary,spkt_elec,condition_elec,[0,0.3],'Elecresponse');
% SPKsummary.spkt_laser=SPKsummary.spkt;
% SPKsummary.spkt_elec=SPKsummary.spkt;
neuron_laser=plot_Spikefiring(SPKsummary,condition_laser,subjectnum,SPKsummary.spkt_laser,[-0.5,1],'firingRate',[0.1,inf],'Layer',Layers(1));
%savefig(['/media/huli/My Book/YLP/TCAdecomposition/Figures/neurons_laser_zscore_amongdepth.fig']);
close gcf;
neuron_laser_correct=cellfun(@(x) basecorrect(x,spkt,-1,-0.2,'subtract'),neuron_laser{1},'UniformOutput',0);
[~,sortindex]=cellfun(@(x) sort(mean(x(spkt>0,:),1)),neuron_laser_correct,'UniformOutput',0);
figure;
for i=1:4
    subplot(1,4,i);
    imagesc(spkt,1:239,neuron_laser_correct{i}(:,sortindex{1})'); axis xy; caxis([-10,20]);xlim([-0.2,0.5]); title(condition_laser{i});
end
colormap jet;
savefig(['/media/huli/My Book/YLP/TCAdecomposition/Figures/neurons_laser_generalplot.fig']);
close gcf;
neuron_elec=plot_Spikefiring(SPKsummary,condition_elec,subjectnum,SPKsummary.spkt_elec,[-0.2,0.5],'firingRate',[0,inf],'Layer',Layers(1));
%savefig(['/media/huli/My Book/YLP/TCAdecomposition/Figures/neurons_elec_zscore_amongdepth.fig']);
close gcf;
neuron_elec_correct=cellfun(@(x) basecorrect(x,spkt,-1,-0.2,'subtract'),neuron_elec{1},'UniformOutput',0);
[~,sortindex]=cellfun(@(x) sort(mean(x(spkt>-0.1&spkt<0.1,:),1)),neuron_elec_correct,'UniformOutput',0);
figure;
for i=1:4
    subplot(1,4,i);
    imagesc(spkt,1:239,neuron_elec_correct{i}(:,sortindex{1})'); axis xy; caxis([-10,20]);xlim([-0.2,0.5]); title(condition_elec{i});
end
colormap jet;
savefig(['/media/huli/My Book/YLP/TCAdecomposition/Figures/neurons_elec_generalplot.fig']);
close gcf;
%% %%%%%%%%%%% deprecated in the manuscript.
[ptmp,ftmp]=timevarying_anova(neuron_laser,SPKsummary.spkt_laser,[-0.5,1]); legend({'P_intensity','P_side','P_interaction'});
eval(['datastat.laser_Pvalue_amongdepth=ptmp;']);
eval(['datastat.laser_Fvalue_amongdepth=ftmp;']);
savefig(['/media/huli/My Book/YLP/TCAdecomposition/Figures/neurons_laser_P_value_correct.fig']);
close gcf;
[ptmp,ftmp]=timevarying_anova(neuron_elec,SPKsummary.spkt_elec,[-0.2,0.5]); legend({'P_intensity','P_side','P_interaction'});
eval(['datastat.elec_Pvalue_amongdepth=ptmp;']);
eval(['datastat.elec_Fvalue_amongdepth=ftmp;']);
savefig(['/media/huli/My Book/YLP/TCAdecomposition/Figures/neurons_elec_P_value_correct.fig']);
close gcf;

% create F value distribution from FiringrateStat.mat
close all;
datastat=matfile('/media/huli/My Book/YLP/TCAdecomposition/Figures/FiringrateStat2.mat');
Savelayer={'All','LayerII_III','LayerIV','LayerVa','LayerVb','LayerVI'};
clear F_* P_*
tmpF_laser=datastat.laser_Fvalue_amongdepth;
tmpF_elec=datastat.elec_Fvalue_amongdepth;
tmpP_laser=datastat.laser_Pvalue_amongdepth;
tmpP_elec=datastat.elec_Pvalue_amongdepth;
for i=1:length(Layers)
   F_laser_intensity(:,i)=tmpF_laser{i}{1};
   F_laser_side(:,i)=tmpF_laser{i}{2};
   F_laser_interaction(:,i)=tmpF_laser{i}{3};
   F_elec_intensity(:,i)=tmpF_elec{i}{1};
   F_elec_side(:,i)=tmpF_elec{i}{2};
   F_elec_interaction(:,i)=tmpF_elec{i}{3};
   P_laser_intensity(:,i)=tmpP_laser{i}{1};
   P_laser_side(:,i)=tmpP_laser{i}{2};
   P_laser_interaction(:,i)=tmpP_laser{i}{3};
   P_elec_intensity(:,i)=tmpP_elec{i}{1};
   P_elec_side(:,i)=tmpP_elec{i}{2};
   P_elec_interaction(:,i)=tmpP_elec{i}{3};
end
figure;
subplot(2,3,1); imagesc(SPKsummary.spkt_laser,1:length(Savelayer),F_laser_intensity'); xlim([-0.5,1]); yticks(1:1:6);yticklabels(Savelayer); axis xy; title('Laser_intensity'); caxis([0,40]);
subplot(2,3,2); imagesc(SPKsummary.spkt_laser,1:length(Savelayer),F_laser_side');xlim([-0.5,1]);yticks(1:1:6);yticklabels(Savelayer); axis xy;title('Laser_side');caxis([0,40]);
subplot(2,3,3); imagesc(SPKsummary.spkt_laser,1:length(Savelayer),F_laser_interaction');xlim([-0.5,1]);yticks(1:1:6);yticklabels(Savelayer); axis xy;title('Laser_interaction');caxis([0,40]);
subplot(2,3,4); imagesc(SPKsummary.spkt_elec,1:length(Savelayer),F_elec_intensity'); xlim([-0.2,0.5]);yticks(1:1:6);yticklabels(Savelayer); axis xy;title('Elec_intensity');caxis([0,40]);
subplot(2,3,5); imagesc(SPKsummary.spkt_elec,1:length(Savelayer),F_elec_side');xlim([-0.2,0.5]);yticks(1:1:6);yticklabels(Savelayer); axis xy;title('Elec_side');caxis([0,40]);
subplot(2,3,6); imagesc(SPKsummary.spkt_elec,1:length(Savelayer),F_elec_interaction');xlim([-0.2,0.5]);yticks(1:1:6);yticklabels(Savelayer); axis xy;title('Elec_interaction');caxis([0,40]);
color=flipud(gray);
colormap(color);
variablenames={'elec','laser'};
factornames={'intensity','side','interaction'};
figure;c=1;
for i=1:length(variablenames)
    for j=1:length(factornames)
        tmp=eval(['P_',variablenames{i},'_',factornames{j}]);
        tmpcorrect=fdr_correct(tmp);
        if i==1; tmpspkt=SPKsummary.spkt_elec; xlimit=[-0.2,0.5];else tmpspkt=SPKsummary.spkt_laser;xlimit=[-0.5,1];end
        subplot(2,3,c); imagesc(tmpspkt,1:length(Savelayer),tmpcorrect');xlim(xlimit);yticks(1:1:6);yticklabels(Savelayer);axis xy; title([variablenames{i},' ',factornames{j}]); caxis([0,0.01]);
        c=c+1;
    end
end%
%
p_laser=cellfun(@(x) timevarying_anova(x,spkt,[-2,4]),neuron_laser,'UniformOutput',0);
p_elec=cellfun(@(x) timevarying_anova(x,spkt,[-2,4]),neuron_elec,'UniformOutput',0);


%
fdr_correct(p_laser,spkt,[-2,4],1,'Laser');
fdr_correct(p_elec,spkt,[-2,4],1,'Elec');
%
fdr_correct(datamat.p_pyramidal_laser,datamat.f_pyramidal_laser,spkt,[-2,4],1,'Pyramidal in Laser');
fdr_correct(datamat.p_interneuron_laser,datamat.f_interneuron_laser,spkt,[-2,4],1,'Interneuron in Laser');
fdr_correct(datamat.p_pyramidal_elec,datamat.f_pyramidal_elec,spkt,[-2,4],1,'Pyramidal in Elec');
fdr_correct(datamat.p_interneuron_elec,datamat.f_interneuron_elec,spkt,[-2,4],1,'Interneuron in Elec');
function NormalizedSPK=normalized(SPKsummary,spkt,pcarange,condition,valid)
% modified for neuraltrajectory
% refer to doi.org/10.1016/j.neuron.2013.11.003
% divided FR for each neuron by the maximun variance across conditions.
tmpdata_all=cell(length(valid),1);
for i=1:length(condition)
    tmpdata=eval(['SPKsummary.',condition{i}(1:end-4)]);
    tmpdata=cellfun(@(x) basecorrect(x,spkt,-0.2,-0.1,'subtract'),tmpdata,'UniformOutput',0);
    tmpdata=cellfun(@(x) mean(x(spkt>pcarange(1)&spkt<pcarange(2),:),2),tmpdata,'UniformOutput',0);
    tmpdata_all=cellfun(@(x,y) cat(1,x,y),tmpdata_all,tmpdata,'UniformOutput',0);
end
NormalizedSPK=cellfun(@(x) (x-mean(x))/sqrt(max(x)),tmpdata_all(valid),'UniformOutput',0); 
% softstd normalized will increased the explaination of pca results.
NormalizedSPK=cell2mat(NormalizedSPK');
end
function NormalizedSPK=normalized_shuffle(SPKsummary,spkt,pcarange,condition,valid)
% resample the data with replacement before normalized data
tmpdata_all=cell(length(valid),1);
for i=1:length(condition)
    tmpdata=eval(['SPKsummary.',condition{i}(1:end-4)]);
    tmpdata=cellfun(@(x) datasample(x,size(x,2),2,'Replace',true),tmpdata,'UniformOutput',0);
    tmpdata=cellfun(@(x) basecorrect(x,spkt,-0.2,-0.1,'subtract'),tmpdata,'UniformOutput',0);
    tmpdata=cellfun(@(x) mean(x(spkt>pcarange(1)&spkt<pcarange(2),:),2),tmpdata,'UniformOutput',0);
    tmpdata_all=cellfun(@(x,y) cat(1,x,y),tmpdata_all,tmpdata,'UniformOutput',0);
end
NormalizedSPK=cellfun(@(x) (x-mean(x))/sqrt(max(x)),tmpdata_all(valid),'UniformOutput',0); 
% softstd normalized will increased the explaination of pca results.
NormalizedSPK=cell2mat(NormalizedSPK');
invalid=isnan(NormalizedSPK(1,:));
NormalizedSPK(:,invalid)=[];
end
function dist=pca_shuffle_baseline(SPKsummary,spkt,pcarange,condition,valid,compare,npcs)
% resample the data with replacement across different conditions
tmpdata_all=cell(length(valid),1);
for c=1:length(compare)
    conditiontmp=condition(compare{c});conditionlength=[];
for i=1:length(condition(compare{c}))
    tmpdata=eval(['SPKsummary.',condition{compare{c}(i)}(1:end-4)]);
    tmpdata_all=cellfun(@(x,y) cat(2,x,y),tmpdata_all,tmpdata,'UniformOutput',0);
    conditionlength{i}=cellfun(@(x) size(x,2), tmpdata,'UniformOutput',1);
end
conditionlength1=sum(cell2mat(conditionlength),2);
labels=[];
for i=1:length(conditionlength1)
    labelstmp=randperm(conditionlength1(i));
    for j=1:length(compare{c})
        labels{i,j}=labelstmp(1:conditionlength{j}(i));
        labelstmp(1:conditionlength{j}(i))=[];
    end
end
for i=1:length(compare{c})
    labeltmp=labels(:,i);
    tmpdata_all2(:,i)=cellfun(@(x,y) x(:,y),tmpdata_all,labeltmp,'UniformOutput',0);
end
tmpshuffledata_all=cell(length(valid),1);
tmpdata=cellfun(@(x,y) datasample(x,length(y),2,'Replace',true),tmpdata_all2,labels,'UniformOutput',0);
tmpdata=cellfun(@(x) mean(x(spkt>pcarange(1)&spkt<pcarange(2),:),2),tmpdata,'UniformOutput',0);
for i=1:length(condition(compare{c}))
  tmpshuffledata_all=cellfun(@(x,y) cat(1,x,y),tmpshuffledata_all,tmpdata(:,i),'UniformOutput',0);
end
NormalizedSPK=cellfun(@(x) (x-mean(x))/sqrt(max(x)),tmpshuffledata_all(valid),'UniformOutput',0); 
% softstd normalized will increased the explaination of pca results.
NormalizedSPK=cell2mat(NormalizedSPK');
invalid=isnan(NormalizedSPK(1,:));
NormalizedSPK(:,invalid)=[];
[~,score,latent,~,explain]=pca(NormalizedSPK);
[score,chopSPKsummary]=chopPCAmatrix(score,NormalizedSPK,sum(spkt>pcarange(1)&spkt<pcarange(2)));
npcs=min(find(cumsum(explain)>85));
[dist(:,c)]=cal_distangle(score,npcs,{1:length(compare{c})});
end
end

function SPKsummary=Rearrange_trialorder(SPKsummary,Highcondition,Lowcondition)
% rearrange the trialorder of the conditions according to the behavior
% score
highscore=eval(['SPKsummary.',Highcondition,'_behaviorscore']);
lowscore=eval(['SPKsummary.',Lowcondition,'_behaviorscore']);
highdata=eval(['SPKsummary.',Highcondition]);
lowdata=eval(['SPKsummary.',Lowcondition]);
hightrialnum=cellfun(@(x) size(x,2),highdata,'UniformOutput',0);
lowtrialnum=cellfun(@(x) size(x,2),lowdata,'UniformOutput',0);
Data=cellfun(@(x,y) cat(2,x,y),highdata,lowdata,'UniformOutput',0);
Datascore=cellfun(@(x,y) cat(1,x,y),highscore,lowscore,'UniformOutput',0);
[~,index]=cellfun(@(x) sort(x),Datascore,'UniformOutput',0);
valid=cellfun(@(x) ~isempty(x),index,'UniformOutput',1);
invalid=cellfun(@(x) isempty(x),index,'UniformOutput',1); % rat6 lack of behaviorscores
newLowdata=cellfun(@(x,y,z) x(:,z(1:y)),Data(valid),lowtrialnum(valid),index(valid),'UniformOutput',0);
newHighdata=cellfun(@(x,y,z) x(:,z(end-y:end)),Data(valid),hightrialnum(valid),index(valid),'UniformOutput',0);
newLowdata=[lowdata(~valid);newLowdata];
newHighdata=[highdata(~valid);newHighdata];
eval(['SPKsummary.',Highcondition,'=newHighdata;']);
eval(['SPKsummary.',Lowcondition,'=newLowdata;']);
    
end
function output=LoadBehaviorScore(filename,condition,subject,output)
%% load behavior score from the xlsx sheet
[a,b]=xlsread(fullfile(subject,filename));
if ~isempty(output.dataoutput)
    index=ismember(b(:,1),{strrep(condition,'_S1_3.mat',[])});
    score=a(index,3);
    output.behaviorscore=repmat({score},[size(output.dataoutput,1),1]);
else
    output.behaviorscore=[];
end
end
function fdr_correct(pSPKsummary,fSPKsummary,lfpt,timeband,option,figtitle)
yaxislabel={'Layer II/III','Layer IV','Layer Va','Layer Vb','Layer VI'};
if option
    for j=1:size(pSPKsummary,2)
        P_intensitytmp=[];P_sidetmp=[];P_interactiontmp=[];
        for k=1:size(timeband,1)
            P_intensitytmp=cat(1,P_intensitytmp,mafdr(pSPKsummary{j}{1}(lfpt>timeband(k,1)&lfpt<timeband(k,2),1),'BHFDR',true));
            P_sidetmp=cat(1,P_sidetmp,mafdr(pSPKsummary{j}{2}(lfpt>timeband(k,1)&lfpt<timeband(k,2),1),'BHFDR',true));
            P_interactiontmp=cat(1,P_interactiontmp,mafdr(pSPKsummary{j}{3}(lfpt>timeband(k,1)&lfpt<timeband(k,2),1),'BHFDR',true));
        end
        P_intensity(:,j)=P_intensitytmp;
         f_intensity(:,j)=fSPKsummary{j}{1};
        P_side(:,j)=P_sidetmp;
        f_side(:,j)=fSPKsummary{j}{2};
        P_interaction(:,j)=P_interactiontmp;
         f_interaction(:,j)=fSPKsummary{j}{3};
    end
    index=zeros(1,length(lfpt));
    for k=1:size(timeband,1)
        index=index|lfpt>timeband(k,1)&lfpt<timeband(k,2);
    end
    lfpt=lfpt(index);
    timeband=reshape(timeband,[],1);
    timeband=[min(timeband),max(timeband)];
    figure; suptitle(figtitle);subplot(2,3,1);imagesc(lfpt,1:size(pSPKsummary,2),P_intensity');xlabel('Time (s)');xlim(timeband);caxis([0,0.05]);set(gca,'YTick',[1:5],'YTickLabel',yaxislabel);title('Intensity P');
    subplot(2,3,2);imagesc(lfpt,1:size(pSPKsummary,2),P_side');xlabel('Time (s)');xlim(timeband);caxis([0,0.05]); set(gca,'YTick',[1:5],'YTickLabel',yaxislabel);title('Side P'); 
    subplot(2,3,3);imagesc(lfpt,1:size(pSPKsummary,2),P_interaction');xlabel('Time (s)');xlim(timeband);caxis([0,0.05]); set(gca,'YTick',[1:5],'YTickLabel',yaxislabel);title('Side*Intensity P'); 
else
    for i=1:size(pSPKsummary,2)
        P_intensity(:,i)=pSPKsummary{i}{1};
        f_intensity(:,i)=fSPKsummary{i}{1};
        P_side(:,i)=pSPKsummary{i}{2};
        f_side(:,i)=fSPKsummary{i}{2};
        P_interaction(:,i)=pSPKsummary{i}{3};
        f_interaction(:,i)=fSPKsummary{i}{3};
    end
    timeband=reshape(timeband,[],1);
    timeband=[min(timeband),max(timeband)];
    figure; suptitle(figtitle);subplot(2,3,1);imagesc(lfpt,1:size(pSPKsummary,2),P_intensity');xlabel('Time (s)');xlim(timeband);caxis([0,0.05]); set(gca,'YTick',[1:5],'YTickLabel',yaxislabel);title('Intensity P'); 
    subplot(2,3,2);imagesc(lfpt,1:size(pSPKsummary,2),P_side');xlabel('Time (s)');xlim(timeband);caxis([0,0.05]); set(gca,'YTick',[1:5],'YTickLabel',yaxislabel);title('Side P'); 
    subplot(2,3,3);imagesc(lfpt,1:size(pSPKsummary,2),P_interaction');xlabel('Time (s)');xlim(timeband);caxis([0,0.05]); set(gca,'YTick',[1:5],'YTickLabel',yaxislabel);title('Side*Intensity P');
end
    subplot(2,3,4);imagesc(lfpt,1:size(pSPKsummary,2),f_intensity');xlim(timeband);caxis([0,20]);title('Intensity P'); set(gca,'YTick',[1:5],'YTickLabel',yaxislabel);set(gca,'FontSize',14);
    subplot(2,3,5);imagesc(lfpt,1:size(pSPKsummary,2),f_side');xlim(timeband);caxis([0,20]);title('Side P'); set(gca,'YTick',[1:5],'YTickLabel',yaxislabel);set(gca,'FontSize',14);
    subplot(2,3,6);imagesc(lfpt,1:size(pSPKsummary,2),f_interaction');xlim(timeband);caxis([0,20]);title('Side*Intensity P'); set(gca,'YTick',[1:5],'YTickLabel',yaxislabel);set(gca,'FontSize',14);
end
function output=Loaddata(path,file,name,propname,prop)
dataoutput=[];spikename=[];spikesubject=[];
try 
    tmp=load(fullfile(path,file),name);
    tmp=eval(['tmp.',name]);
    data=getfield(tmp,'rawdata');
    for i=1:length(data)
        for j=1:length(data{i})
        data{i}{j}(data{i}{j}>-0.002&data{i}{j}<0.002)=[]; % delete electrical induced noise around 0.
        end
    end
     chooseinfo=getfield(tmp,'Chooseinfo');% get the spikename
switch propname
    case 'rawtobin'
        for i=1:length(data)
        [binoutput,t]=cellfun(@(x) binspikes(x,1/prop,[-2,4]),data{i},'UniformOutput',0);
        dataoutput=cat(1,dataoutput,{mean(cell2mat(binoutput),2)});
        spikename=cat(1,spikename,chooseinfo.spikename(i));
        spikesubject=cat(1,spikesubject,{name});
        output.spkt=t{1};
        end   
    case 'rawtobin-noaverage'
        for i=1:length(data)
        [binoutput,t]=cellfun(@(x) binspikes(x,1/prop,[-2,4]),data{i},'UniformOutput',0);
        dataoutput=cat(1,dataoutput,{cell2mat(binoutput)});
        spikename=cat(1,spikename,chooseinfo.spikename(i));
        spikesubject=cat(1,spikesubject,{name});
        output.spkt=t{1};
        end   
    case 'rawtopsth'
        for i=1:length(data)
        [binoutput,spkt]=cal_psth(data{i},prop,'n',[-2,4],1);
            if isempty(binoutput)
            binoutput=nan(1,length(spkt));
            end
            binoutput=binoutput';
            output.spkt=spkt;
        dataoutput=cat(1,dataoutput,{binoutput});
        spikename=cat(1,spikename,chooseinfo.spikename(i));
        spikesubject=cat(1,spikesubject,{name});
        end   
    case 'rawtopsth-noaverage'       
        for i=1:length(data)
            binoutputall=[];binoutput=[];
            for j=1:size(data{i},2)      
        [binoutputtmp,spkt]=cal_psth(data{i}(j),prop,'n',[-2,4],1);
            if isempty(binoutputtmp)
            binoutputtmp=zeros(1,length(spkt));
            end
            binoutput(:,j)=binoutputtmp';
            end
            binoutputall=binoutput;
            output.spkt=spkt;
        dataoutput=cat(1,dataoutput,{binoutputall});
        spikename=cat(1,spikename,chooseinfo.spikename(i));
        spikesubject=cat(1,spikesubject,{name});
        end   
    case 'binneddata'
            for i=1:length(data)
                dataoutput=cat(1,dataoutput,{mean(data{i},2)});
                spikename=cat(1,spikename,chooseinfo.spikename(i));
                spikesubject=cat(1,spikesubject,{name});
            end
    case 'rawraster'
          for i=1:length(data)
                [binoutput{i}]=cell2mat(cellfun(@(x) binspikes(x,1000,[-2,4]),data{i},'UniformOutput',0));
          end
        dataoutput=cat(2,dataoutput,cellfun(@(x) sparse(x),binoutput,'UniformOutput',0));
        output.spkt=linspace(-2,4,6001);
        spikename=cat(1,spikename,chooseinfo.spikename);
        spikesubject=cat(1,spikesubject,repmat({name},[size(dataoutput,2),1])); 
         dataoutput=dataoutput';
      end
catch
            dataoutput=cat(1,dataoutput,{});
            spikename=cat(1,spikename,{});
            spikesubject=cat(1,spikesubject,{});
         end
output.dataoutput=dataoutput;
output.spikename=spikename;
output.spikesubject=spikesubject;

end
function output=Loaddata_suffle(path,file_all,name,propname,prop,shuffle)
% used for pca suffle data from elec and laser modality
data_all=[];con=1;eventindex=[];n=1;
 for v=1:length(file_all)
    for c=1:length(file_all{v})
        tmp=load(fullfile(path{v},file_all{v}{c}),name);
        tmp=eval(['tmp.',name]);
        data=getfield(tmp,'rawdata');
        for i=1:length(data) %spike
            if v==1&&c==1
                data_all{i}=[];
            end
            for j=1:length(data{i})%trial
            data{i}{j}(data{i}{j}>-0.002&data{i}{j}<0.002)=[]; % delete electrical induced noise around 0.
            end
            data_all{i}=cat(2,data_all{i},data{i});
        end
        chooseinfo=getfield(tmp,'Chooseinfo');% get the spikename
        conditionlength(con)=length(chooseinfo.Eventindex);
        eventindex=cat(1,eventindex,repmat(file_all{v}(c),[conditionlength(con),1]));
        conditionname{con}=file_all{v}{c};
        con=con+1;
    end
 end
for m=1:shuffle
%     indices=cvpartition(eventindex,'KFold',length(conditionlength),'Stratify',true);
%     for i=1:length(conditionlength)
%         eventindexresample{i}=false(length(eventindex),1);
%         eventindexresample{i}=indices.test(i);
%     end
%      indices=datasample(eventindex,length(eventindex),'Replace',false);
%      for i=1:length(conditionlength)
%          eventindexresample{i}=strcmp(indices,conditionname{i});
%      end
  %indices=datasample(eventindex,length(eventindex),'Replace',true);
for c=1:length(conditionname)
    index=ismember(eventindex,conditionname{c});
    data=cellfun(@(x) x(index),data_all,'UniformOutput',0);
    data=cellfun(@(x) datasample(x,length(x),'Replace',true),data,'UniformOutput',0);
%     data=cellfun(@(x) x(eventindexresample{c}),data_all_tmp,'UniformOutput',0);
    dataoutput=[]; spikesubject=[];spikename=[];
switch propname
    case 'rawtobin'
        for i=1:length(data)
        [binoutput,t]=cellfun(@(x) binspikes(x,1/prop,[-2,4]),data{i},'UniformOutput',0);
        dataoutput=cat(1,dataoutput,{mean(cell2mat(binoutput),2)});
        spikename=cat(1,spikename,chooseinfo.spikename(i));
        spikesubject=cat(1,spikesubject,{name});
        output(m).spkt=t{1};
        end   
    case 'rawtobin-noaverage'
        for i=1:length(data)
        [binoutput,t]=cellfun(@(x) binspikes(x,1/prop,[-2,4]),data{i},'UniformOutput',0);
        dataoutput=cat(1,dataoutput,{cell2mat(binoutput)});
        spikename=cat(1,spikename,chooseinfo.spikename(i));
        spikesubject=cat(1,spikesubject,{name});
        output(m).spkt=t{1};
        end   
    case 'rawtopsth'  
        for i=1:length(data)
            [binoutput,spkt]=cal_psth(data{i},prop,'n',[-2,4],1);
            if isempty(binoutput)
            binoutput=nan(1,length(spkt));
            end
            binoutput=binoutput';
            output(m).spkt=spkt;
        dataoutput=cat(1,dataoutput,{binoutput});
        spikename=cat(1,spikename,chooseinfo.spikename(i));
        spikesubject=cat(1,spikesubject,{name});
        end   
    case 'rawtopsth-noaverage'      
        for i=1:length(data_all)
            binoutputall=[];binoutput=[];
            for j=1:size(data{i},2)      
        [binoutputtmp,spkt]=cal_psth(data{i}(j),prop,'n',[-2,4],1);
            if isempty(binoutputtmp)
            binoutputtmp=zeros(1,length(spkt));
            end
            binoutput(:,j)=binoutputtmp';
            end
            binoutputall=binoutput;
            output(m).spkt=spkt;
        dataoutput=cat(1,dataoutput,{binoutputall});
        spikename=cat(1,spikename,chooseinfo.spikename(i));
        spikesubject=cat(1,spikesubject,{name});
        end   
    case 'binneddata'
            for i=1:length(data)
                dataoutput=cat(1,dataoutput,{mean(data{i},2)});
                spikename=cat(1,spikename,chooseinfo.spikename(i));
                spikesubject=cat(1,spikesubject,{name});
            end
    case 'rawraster'
          for i=1:length(data)
                [binoutput{i}]=cell2mat(cellfun(@(x) binspikes(x,40000,[-2,4]),data{i},'UniformOutput',0));
          end
        dataoutput=cat(2,dataoutput,cellfun(@(x) sparse(x),binoutput,'UniformOutput',0));
        output(m).spkt=linspace(-2,4,240001);
        spikename=cat(1,spikename,chooseinfo.spikename(i));
        spikesubject=cat(1,spikesubject,{name}); 
        dataoutput=dataoutput';
end
    output(m).spikename=spikename;
    output(m).spikesubject=spikesubject;
    eval(['output(m).',conditionname{c}(1:end-4),'=dataoutput;']);
    end
 end
end
function [outputdata,t]=cal_psth(data,sig,plt,T,err)
for i=1:length(data)
    output(i).times=data{i};
end
[outputdata,t]=psth(output,sig,plt,T,err);
end
function output=catall_shuffle(shuffle_data,conditionname,shuffleindex)

    for c=1:length(shuffle_data)
         output.spikename=[];
         output.spikesubject=[];
        for i=1:length(conditionname)
            if c==1;eval(['output.',conditionname{i}(1:end-4),'=[];']);end
            eval(['output.',conditionname{i}(1:end-4),'=cat(1,output.',conditionname{i}(1:end-4),',shuffle_data{c}(shuffleindex).',conditionname{i}(1:end-4),');']);
        end
        output.spikename=cat(1,output.spikename,shuffle_data{c}.spikename);
        output.spikesubject=cat(1,output.spikesubject,shuffle_data{c}.spikesubject);
        end
end
function data=Relabel_shuffle(data,conditionname)
    
    spikedata=[];
    for i=1:length(conditionname)
        tmp=eval(['data.',conditionname{i}(1:end-4),';']);
        if i==1
            neuronnumber=length(tmp);
        end
        spikedata=cat(2,spikedata,tmp);
    end
    for i=1:neuronnumber
        for c=1:length(conditionname)
            [index]=randperm(4);
            shuffledata(i,c)=spikedata(i,index(c));
        end
    end
    for i=1:length(conditionname)
        eval(['data.',conditionname{i}(1:end-4),'=shuffledata(:,i);']);
    end
end
function output=catall(data)
output.spikename=[];
output.dataoutput=[];
output.spikesubject=[];
output.behaviorscore=[];
for i=1:length(data)
    if ~isempty(data{i}.dataoutput)
    output.dataoutput=cat(1,output.dataoutput,data{i}.dataoutput);
    output.spikename=cat(1,output.spikename,data{i}.spikename);
    output.spikesubject=cat(1,output.spikesubject,data{i}.spikesubject);
    output.behaviorscore=cat(1,output.behaviorscore,data{i}.behaviorscore);
    end
end
try
    output.spkt=data{2}.spkt;
end
end
function output=spikeproperties(output)
% write properties of each spike from the cellmetric.mat
% including firingrate, thoughtoPeak,abratio,putativeCelltype
properties={'firingRate','troughToPeak','putativeCellType','ab_ratio','maxWaveformCh1'};
for i=1:length(output.spikesubject)
    try
    cellmetrics=matfile([output.spikesubject{i},'/',output.spikesubject{i},'.cell_metrics.cellinfo.mat']);
    cellmetrics=cellmetrics.cell_metrics;
    clu=str2num(strrep(output.spikename{i},'cluster1_',[]));
    index=find(cellmetrics.cluID==clu);
    for j=1:length(properties)
    eval(['output.',properties{j},'(i)=cellmetrics.',properties{j},'(index);']);
    end
    catch
    for j=1:length(properties)
        eval(['output.',properties{j},'(i)=[];']);
    end 
    end
end
end
function len=channellength(x)
if ischar(x)
    len=length(str2num(x));
elseif isnumeric(x)
    if isnan(x)
        len=0;
    else
        len=length(x);
    end
end
end
function adddepthline(timerange)
  separate=[1.5,5.5,8.5,10.5,13.5];
    hold on;line(gca,[timerange(1),timerange(2)],[separate(1),separate(1)],'LineWidth',1); % separate layer I and II/III
    line(gca,[timerange(1),timerange(2)],[separate(2),separate(2)],'LineWidth',1); % separate layer II/III and layer IV
      line(gca,[timerange(1),timerange(2)],[separate(3),separate(3)],'LineWidth',1); % separate layer IV and layer Va
       line(gca,[timerange(1),timerange(2)],[separate(4),separate(4)],'LineWidth',1,'LineStyle','--'); % separate layer Va and layer Vb
       line(gca,[timerange(1),timerange(2)],[separate(5),separate(5)],'LineWidth',1); % separate layer Vb and layer VI
       set(gca,'YTick',[separate(1)/2,separate(1)+(separate(2)-separate(1))/2, separate(2)+(separate(3)-separate(2))/2,separate(3)+(separate(4)-separate(3))/2,separate(4)+(separate(5)-separate(4))/2,separate(5)+(16-separate(5))/2]);
       set(gca,'YTickLabel',{'Layer I','Layer II/III','Layer IV','Layer Va','Layer Vb','Layer VI'});
       xlabel('troughToPeak');
end
function index=conditionfilter(data,varargin)
p=inputParser;
addParameter(p,'putativeCellType',[]);
addParameter(p,'Layer',[]);
addParameter(p,'firingRate',[]);
addParameter(p,'invalid',[]);
parse(p,varargin{:}{:});
indexType=true(1,length(data.Layer));
indexLayer=true(1,length(data.Layer));
indexfiring=true(1,length(data.Layer));
indexinvalid=true(1,length(data.Layer));
binplotoutput=[];
if ~isempty(p.Results.putativeCellType)
    indexType=ismember(data.putativeCellType,p.Results.putativeCellType);
end
if ~isempty(p.Results.Layer)
    for j=1:length(p.Results.Layer)
        if ~strcmp(p.Results.Layer{j},'All')
        indexLayer(j,:)=ismember(data.Layer,p.Results.Layer{j});
        else
            indexLayer(j,:)=true(size(data.Layer));
        end
    end
end
if ~isempty(p.Results.firingRate)
    indexfiring=data.firingRate>p.Results.firingRate(1)&data.firingRate<p.Results.firingRate(2);
end
if ~isempty(p.Results.invalid)
    indexinvalid=~p.Results.invalid;
end
if ~isempty(p.Results.Layer)
for j=1:length(p.Results.Layer)
    index{j}=indexType&indexfiring&indexinvalid&indexLayer(j,:);
end
else
    index{1}=indexType&indexfiring&indexinvalid;
end
end
function binplotoutput=plot_Spikefiring(data,condition,subjectnum,spkt,xrange,varargin)
p=inputParser;
addParameter(p,'putativeCellType',[]);
addParameter(p,'Layer',[]);
addParameter(p,'firingRate',[]);
addParameter(p,'invalid',[]);
parse(p,varargin{:});
index=conditionfilter(data,varargin);
figure; conditioncolor={[0,0.447,0.741],[0.85,0.325,0.098],[0.929,0.694,0.125],[0.494,0.184,0.556],'black'};
  if isempty(p.Results.Layer)
      Layer={'All'};
  else
      Layer=p.Results.Layer;
  end
 for j=1:length(Layer) 
    subplot(1,length(index),j);
    invalidindex=[];
    for i=1:length(condition)
    bindata=eval(['data.',condition{i}(1:end-4)]);
    bin_plot=bindata(index{j});
    binplot=[];
    for k=1:length(bin_plot)
        if ~isempty(bin_plot{k})
%             if ~isnan(bin_plot{k})
            binplottmp=basecorrect(nanmean(bin_plot{k},2),spkt,-0.5,0,'zscore');
            binplot=cat(2,binplot,binplottmp);
%             end
        end
    end
   
%     shadebar(spkt,binplot,conditioncolor{i});
    
    binplotoutput{j}{i}=cell2mat(cellfun(@(x) mean(x,2),bin_plot,'UniformOutput',0)');
    invalid=isnan(binplot(1,:))|isinf(binplot(1,:));
    tmpinvalid=1:size(binplot,2);
    tmpinvalid(invalid)=[];
    binplot=binplot(:,~invalid);
    [~,invalid2]=deleteoutliers(max(binplot));
    binplot(:,invalid2)=[];
    tmpinvalid(invalid2)=[];
    binplot_plot=nanmean(binplot,2);
     disp(['Number of',Layer{j},'_',p.Results.putativeCellType,' is:',num2str(size(binplot,2))]);
    %plot(spkt,smooth(binplot_plot,3),'color',conditioncolor{i},'LineWidth',1.5);
    plot(spkt,smooth(binplot_plot,3),'color',conditioncolor{i},'LineWidth',1.5);
%     shadebar(spkt,binplot,conditioncolor{i});
    hold on;xlim(xrange);xlabel('Time (s)'); ylabel('Firing Rate(Hz)'); legend(condition);
    invalidindex=cat(2,invalidindex,tmpinvalid);
    end
    tmpinvalid=tabulate(invalidindex);
    validindex=tmpinvalid(tmpinvalid(:,2)==length(condition),1);
    binplotoutput{j}=cellfun(@(x) x(:,validindex),binplotoutput{j},'UniformOutput',0);
    title([Layer{j},'N=',num2str(length(validindex))]);
end

end
function SPKsummary=cal_responses(SPKsummary,spkt,conditionname,Roit,option)
groupall={'Intensity','Side'};
groupcom={[1,1],[1,2],[2,1],[2,2]};
for i=1:length(conditionname)
    tmpmat=eval(['SPKsummary.',conditionname{i}(1:end-4)]);
    for j=1:length(tmpmat)
        for k=1:size(tmpmat{j},2)
            tmpsvd{j}{i}(k)=mean(tmpmat{j}(spkt>Roit(1)&spkt<Roit(2),k),1);
            tmpgroup{j}{i}(:,k)=groupcom{i};
        end
    end
end
for j=1:length(tmpsvd)
    tmpstatsvd{j}=cell2mat(tmpsvd{j});
    tmpstatgroup{j}=cell2mat(tmpgroup{j});
    [~,tbl]=anovan(tmpstatsvd{j},{tmpstatgroup{j}(1,:),tmpstatgroup{j}(2,:)},'model','interaction','varnames',{'Intensity','Side'},'display','off');
    tmpstr='';
    if tbl{2,7}<0.05 
        tmpstr=[tmpstr,'Intensity'];
    end
    if tbl{3,7}<0.05
        tmpstr=[tmpstr,'Side'];
    end   
    if isempty(tmpstr)&&tbl{4,7}<0.05
       tmpstr='Interaction';  
    elseif isempty(tmpstr) 
        tmpstr='none';
    end
    if strcmp(tmpstr,'IntensitySide')
        tmpstr='Both';
    end
    eval(['SPKsummary.',option,'{j}=tmpstr;']);
end
        
end
function [pSPK,fSPK]=timevarying_anova(SPKsummarydata_all,spkt,timerange)
%% calculate the two way anova among depth of SPK data intensity&side 2*2 in each time point
for c=1:length(SPKsummarydata_all)
    SPKsummarydata=SPKsummarydata_all{c}; datamat=[];
F_intensity=nan(size(SPKsummarydata{1},1),1);
P_intensity=nan(size(SPKsummarydata{1},1),1);
F_Side=nan(size(SPKsummarydata{1},1),1);
P_Side=nan(size(SPKsummarydata{1},1),1);
F_interaction=nan(size(SPKsummarydata{1},1),1);
P_interaction=nan(size(SPKsummarydata{1},1),1);
for i=1:size(SPKsummarydata{1})
    nanindex=isnan(SPKsummarydata{1}(1,:))|isnan(SPKsummarydata{2}(1,:))|isnan(SPKsummarydata{3}(1,:))|isnan(SPKsummarydata{4}(1,:));
    if spkt(i)>=timerange(1)&&spkt(i)<=timerange(2)
        datamat(:,1,1)=SPKsummarydata{1}(i,~nanindex);
        datamat(:,1,2)=SPKsummarydata{2}(i,~nanindex);
        datamat(:,2,1)=SPKsummarydata{3}(i,~nanindex);
        datamat(:,2,2)=SPKsummarydata{4}(i,~nanindex);
        try
            [tbl,rm] = simple_mixed_anova(datamat, [],{'Intensity','Side'});
            F_intensity(i)=tbl{3,4};
            P_intensity(i)=tbl{3,5};
            F_Side(i)=tbl{5,4};
            P_Side(i)=tbl{5,5};
            F_interaction(i)=tbl{7,4};
            P_interaction(i)=tbl{7,5}; 
        catch
            F_intensity(i)=NaN;
            P_intensity(i)=NaN;
            F_Side(i)=NaN;
            P_Side(i)=NaN;
            F_interaction(i)=NaN;
            P_interaction(i)=NaN; 
        end
    end
end
% figure; subplot(1,3,1);imagesc(lfpt,1:16,P_intensity');xlim([0,0.1]);caxis([0,0.05]);title('Intensity P');
% subplot(1,3,2);imagesc(lfpt,1:16,P_Side');xlim([0,0.1]);caxis([0,0.05]);title('Side P');
% subplot(1,3,3);imagesc(lfpt,1:16,P_interaction');xlim([0,0.1]);caxis([0,0.05]);title('Side*Intensity P');
 pSPK{c}={P_intensity,P_Side,P_interaction};
 fSPK{c}={F_intensity,F_Side,F_interaction};
 subplot(1,length(SPKsummarydata_all),c);
 P_intensity_correct=mafdr(P_intensity(spkt>=timerange(1)&spkt<=timerange(2)),'BHFDR',true);
 P_side_correct=mafdr(P_Side(spkt>=timerange(1)&spkt<=timerange(2)),'BHFDR',true);
 P_interaction_correct=mafdr(P_interaction(spkt>=timerange(1)&spkt<=timerange(2)),'BHFDR',true);
 plot(spkt(spkt>=timerange(1)&spkt<=timerange(2)),P_intensity_correct,'r');
 hold on;
 plot(spkt(spkt>=timerange(1)&spkt<=timerange(2)),P_side_correct,'b');
 plot(spkt(spkt>=timerange(1)&spkt<=timerange(2)),P_interaction_correct,'g');
 xlim([-0.2,0.5]); ylim([0,0.05]);
end
end

function output=tabulate_all(data,categroy)
% tabulate the data according to the given categroy
output=cell(length(categroy),3);
output(:,1)=categroy;
output(:,2)=num2cell(zeros(length(categroy),1));
output(:,3)=num2cell(zeros(length(categroy),1));
outputtmp=tabulate(data);
for i=1:length(categroy)
    if sum(ismember(outputtmp(:,1),categroy(i)))==1
        output(i,:)=outputtmp(ismember(outputtmp(:,1),categroy(i)),:);
    end
end
end
function spikepresent_show(SPKsummary, spkt, condition, index,timelimit)
    for i=1:1:length(condition)
        for c=1:length(index)
        tmpspk=eval(['SPKsummary.',condition{i}(1:end-4),'{index(c)};']);  
%         if i==3; tmpspk(:,6)=[]; end;if i==1; tmpspk(:,3)=[];end
        if ~issparse(tmpspk) 
      
%         tmpspk_all(:,c)=basecorrect(tmpspk,spkt,-0.1,0,'Subtract');
%         tmpspk=basecorrect(mean(tmpspk,2),spkt,-0.1,0,'Zscore');
        tmpspk_all(:,c)=smoothdata(mean(tmpspk,2),'gaussian',5);  
        tmpspk_all(:,c)=basecorrect(tmpspk_all(:,c),spkt,-2,0,'Zscore');  
        hold on; plot(spkt,mean(tmpspk_all,2)); xlim(timelimit);
        elseif c==1
            subplot(2,2,i)
            [~,xPoints,yPoints]=plotSpikeRaster(logical(tmpspk)','PlotType','vertline2','TimePerBin',1/40000);
            plot(gca,xPoints/40000-2,yPoints);  axis  tight; xlim(timelimit);
            title(condition{i});
        end
        end
              
    end
end
function [coeff,SPKdata]=chopPCAmatrix(score,SPKsummarydata,choplength)
choplength=repmat(choplength,[size(score,1)/choplength,1]);
    for i=1:length(choplength)
        coeff(:,:,i)=score(1:choplength(i),:);
        SPKdata(:,:,i)=SPKsummarydata(1:choplength(i),:);
        score(1:choplength(i),:)=[];
        SPKsummarydata(1:choplength(i),:)=[];
    end
end
function [dist,degree,ang]=cal_distangle(score,npcs,compare)  
    score=score(:,1:npcs,:);

%    distlength_all=norm(score(:,:,1)-score(:,:,2))+norm(score(:,:,1)-score(:,:,3))+norm(score(:,:,1)-score(:,:,4))+norm(score(:,:,2)-score(:,:,3))+norm(score(:,:,2)-score(:,:,4))+norm(score(:,:,3)-score(:,:,4));
%distlength_all=1;
 for i=1:length(compare)
     if length(compare{i})==2
    comp1=score(:,:,compare{i}(1));
    comp2=score(:,:,compare{i}(2));
    %dist(i)=norm(comp1-comp2)/norm(comp1+comp2);
%     for j=1:size(comp1,1) dist(j,i)=norm(comp1(j,:)-comp2(j,:));
%     end
    %dist(i)=norm(comp1-comp2)/distlength_all;
    [dist(:,i),angle(:,i)]=meculidean(comp1,comp2);
    %dist(i)=seculidean(comp1,comp2);
    angle(i)=subspace(comp1,comp2);
    procrute(i)=procrustes(comp1,comp2,'reflection',false);
     else
         subcompare={[compare{i}(1),compare{i}(2)],[compare{i}(1),compare{i}(3)]};
         subdist=[];
         for c=1:2
            comp1=score(:,:,subcompare{c}(1));
            comp2=score(:,:,subcompare{c}(2));
            [subdist(:,c),subang(:,c)]=meculidean(comp1,comp2);
         end
         dist(:,i)=(subdist(:,1)-subdist(:,2));
         ang(:,i)=(subang(:,1)-subang(:,2));
         degree(:,i)=meculidean(comp1,comp1(1:10,:));
     end 
 end

 

end
function [d,ang]=meculidean(pt1,pt2)
for i=1:size(pt1,1)
    disttmp=squeeze(sqrt(sum((pt1(i,:)'-pt2').^2)));
    [d(i),index]=min(disttmp);
    ang(i)=acos(dot(pt1(i,:)./norm(pt1(i,:)),pt2(index,:)./norm(pt2(index,:))));
end
end
function h=timevaryingttest(data1,data2)
    for i=1:size(data1,2)
        [~,p(i)]=ttest2(data1(:,i),data2(:,i));
    end
    p=fdr_BH(p,0.05);
    h=p<0.01;
end
        