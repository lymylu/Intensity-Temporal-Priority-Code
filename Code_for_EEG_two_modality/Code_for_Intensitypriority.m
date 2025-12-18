% close all
% clc
%%
%ICAremoval;
%%
clear all;
processtypeall={'nose','Fz','_averageref'};
subject=[2:74,76:109];
Location={'R','L'};
Modality={'L','S','V','A'};
basepath='/mnt/Share/yuelp/Data_Fourmodalities/ICAremovalnew';
outputpath='/mnt/Share/yuelp/Data_Fourmodalities/Intensitypriority';
lfpt=linspace(-1,1.999,3000);
for all=2
    processtype=processtypeall{all};% [],'_hightrial'.....
    datamat=matfile(fullfile(outputpath,['trialaverage_',processtypeall{all},'.mat']),'Writable',true);
    EEG=pop_loadset('filename',strcat('S',num2str(subject(1)),'L_ICA.set'),'filepath',basepath);
    if all==2
        EEG=pop_reref(EEG,'Fz');
    end
    datamat.chanlocs=EEG.chanlocs;
    ERPsubject=nan(3000,length(EEG.chanlocs),4,4,107);Specsubject_mwt=nan(3000,6,length(EEG.chanlocs),4,4,107);Specsubject_stft=nan(3000,6,length(EEG.chanlocs),4,4,107);   
   % time(3000)*(freq,100)*modality(L,S,V,A;4)*condition(HL,HR,LL,LR;4)*subject(107)
    
    parfor i=1:length(subject)
         %[ERPLeftHigh,ERPLeftLow,SpecLeftHigh,SpecLeftLow]=getERPandSpecdata(strcat('S',num2str(subject(i)),'L_ICA.set'),basepath,processtype,true,'mwt',true,true);
         %[ERPRightHigh,ERPRightLow,SpecRightHigh,SpecRightLow]=getERPandSpecdata(strcat('S',num2str(subject(i)),'R_ICA.set'),basepath,processtype,true,'mwt',true,true);
        [ERPLeftHigh,ERPLeftLow]=getERPandSpecdata(strcat('S',num2str(subject(i)),'L_ICA.set'),basepath,processtype,true,'mwt',true,false);
       [ERPRightHigh,ERPRightLow]=getERPandSpecdata(strcat('S',num2str(subject(i)),'R_ICA.set'),basepath,processtype,true,'mwt',true,false);
        ERPsubject(:,:,:,:,i)=cat(4,ERPLeftHigh,ERPRightHigh,ERPLeftLow,ERPRightLow);
        %Specsubject_mwt(:,:,:,:,:,i)=cat(5,SpecLeftHigh,SpecRightHigh,SpecLeftLow,SpecRightLow);
%         [~,~,SpecLeftHigh,SpecLeftLow]=getERPandSpecdata(strcat(num2str(subject(i)),'_LH.set'),basepath,processtype,true,'stft',true);
%         [~,~,SpecRightHigh,SpecRightLow]=getERPandSpecdata(strcat(num2str(subject(i)),'_RH.set'),basepath,processtype,true,'stft',true);
%         Specsubject_stft(:,:,:,:,i)=cat(4,SpecLeftHigh,SpecRightHigh,SpecLeftLow,SpecRightLow);
   end
datamat.ERPsubject=ERPsubject;
%datamat.Specsubject_mwt=Specsubject_mwt;
% datamat.Specsubject_stft=Specsubject_stft;
end
%% general plot for laser and electrical modality Fz version
datamat=matfile('trialaverage_Fz.mat');
ERP_EEG=datamat.ERPsubject; % time*channel*modality*condition(hl,hr,ll,lr)*subject
chanlocs_EEG=datamat.chanlocs;
chanlocs=struct2table(chanlocs_EEG);
chanlocs=chanlocs.labels;
channellabel={'C3','C4'};
c=1;
lfpt=linspace(-1,1.999,3000);
modality={'Laser','Elec'}
ERP_EEG=basecorrect(ERP_EEG,lfpt,-0.5,0,'Subtract');
for i=1:2
    for l=1:2
        ERP=ERP_EEG(:,ismember(chanlocs,channellabel(l)),i,:,:);
        subplot(2,2,c)
        plot(lfpt,squeeze(mean(ERP,5))); title([modality{i},' ',channellabel{l}]); xlim([-0.5,1]); set(gca,'YDir','reverse'); ylim([-10,10]);
    for cond=1:4
        %  if all==1 % N2 and P2
        % tindex={lfpt>0.19&lfpt<0.21,lfpt>0.33&lfpt<0.35};
        % elseif all==2 % N1
            if i==1
                tindex={lfpt>0.16&lfpt<0.18};
            else
                tindex={lfpt>0.08&lfpt<0.10};
            end
        % end     
        for tin=1:length(tindex)
        axes('position',[0.1*cond+(l-1)*0.4,1-i*0.5+tin*0.05,0.1,0.1]);   
        topoplot(squeeze(mean(ERP_EEG(tindex{tin},:,i,cond,:),[1,5])),chanlocs_EEG);caxis([-5,5]);
        end
    end
    c=c+1;
    end
end
 
%% Tensor Decomposition from the EEG data
% prepare the data matrix for tensor decomposition.
tensormatrix=matfile(fullfile(outputpath,'tensormatrix_mwt.mat'),'Writable',true);
%invalid=[6,7,12,29,34,35,36,39,40,43,44,56,64,71,73,84,90];
processtype='';
t=linspace(-1,1.999,3000);
Intensity={'High','Low'};
location={'Right','Left'}; % type order RightHigh,RightLow,LeftHigh,LeftLow
Tensormatrix=zeros(2,2,59,sum(t>=-0.2&t<=0.5),6,length(subject),length(Modality));

for c=1:length(subject)
    datamat=matfile(fullfile(basepath,[num2str(subject(c)),'_subject',processtype,'.mat']));
    for m=1:4
    for i=1:2
       for j=1:2
             tmp_spec=evalvar(datamat,[location{i},Intensity{j},'_mwt_',Modality{m}]);
%              tmp_spec=permute(tmp_spec,[2,1,3]);
%              tmp_spec=basecorrect(tmp_spec,t,-0.4,-0.1,'normalized');
             Tensormatrix(i,j,:,:,:,c,m)=permute(squeeze(averagefreq(permute(tmp_spec(:,t>=-0.2&t<=0.5,:),[2,3,1]))),[3,1,2]);
       end
    end
    end
end

%Tensormatrix(:,:,:,:,:,invalid,:)=[];
%  use all 107*4 as the event trial
c=1;
Tensormatrix_all=[];
for i=1:2
    for j=1:2
        Tensormatrix_all=cat(4,Tensormatrix_all,squeeze(Tensormatrix(i,j,:,:,:,:,:)));
        c=c+1;
    end
end
% channel* time *frequency * [condition (LH,LL,RH,RL)*subject]*modality
tensormatrix.Tensormatrix=Tensormatrix_all;
% tensormatrix.Tensormatrix_two=tensormatrix_freqaverage(:,:,:,:,1:2);
% tensormatrix.Tensormatrix_laser=tensormatrix_freqaverage(:,:,:,:,1);
% tensormatrix.Tensormatrix_elec=tensormatrix_freqaverage(:,:,:,:,2);
%%
function [Highdata,Lowdata,HighScore,LowScore]=getEEGdata(filename,basepath,preprocess,modality)
       EEG = pop_loadset('filename',filename,'filepath',basepath);
       if contains(preprocess,'Fz')
        EEG = pop_reref( EEG, 'Fz');
        EEG = eeg_checkset( EEG );
       elseif contains(preprocess,'average')
           EEG=pop_reref(EEG,[]);
           EEG=eeg_checkset(EEG);
       end
       Event=struct2table(EEG.event);
       eventdescription=Event.type;
       modalityindex=contains(eventdescription,modality);
       score=Event.rating;
       invalid=isnan(score(:,1));
       if strcmp(modality,'L')
            invalid=invalid|score(:,1)<4;
       end
       modalityindex(invalid)=false;
       EEG=pop_select(EEG,'trial',find(modalityindex==1));
       Event=struct2table(EEG.event);
       score=Event.rating;
       score=score(:,1);
        [~,index]=sort(score);
        lowscore=index(1:round(length(index)/2));
        highscore=index((length(index)-round(length(index)/2))+1:end);
        lowindex=lowscore;highindex=highscore;
%         if contains(preprocess,'hightrial')
%             lowindex=highscore(1:round(length(highscore)/2));
%             highindex=highscore(length(highscore)-round(length(highscore)/2)+1:end); 
%         end
        HighScore=score(highindex);LowScore=score(lowindex);
        Highdata=mean(EEG.data(:,:,highindex),3);
        Lowdata=mean(EEG.data(:,:,lowindex),3);
        % Leftlow{i}=EEG.data(:,:,ismember(eventdescription,{'S 11','S 12'}));
% %     Lefthigh{i}=EEG.data(:,:,ismember(eventdescription,{'S 14','S 13'}));
end
function [Highdata,Lowdata]=getSpecdata(filename,basepath,preprocess,modality)
      EEG = pop_loadset('filename',filename,'filepath',basepath);
       if contains(preprocess,'Fz')
        EEG = pop_reref( EEG, 'Fz');
        EEG = eeg_checkset( EEG );
       elseif contains(preprocess,'average')
           EEG=pop_reref(EEG,[]);
           EEG=eeg_checkset(EEG);
       end
       Event=struct2table(EEG.event);
       eventdescription=Event.type;
       modalityindex=contains(eventdescription,modality);
       score=Event.rating;
       invalid=isnan(score(:,1));
       if strcmp(modality,'L')
            invalid=invalid|score(:,1)<4;   
       end
       modalityindex(invalid)=false;
       EEG=pop_select(EEG,'trial',find(modalityindex==1));
       Event=struct2table(EEG.event);
       score=Event.rating;
       score=score(:,1);
        [~,index]=sort(score);
        %    % cal the spectrogram.
           specdata=zeros(size(EEG.data,1),3000,100,size(EEG.data,3));
           for c=1:size(EEG.data,1)
               for d=1:size(EEG.data,3)
%                    if contains(preprocess,'FFT')
                      %[~, tmp] = sub_stft(squeeze(EEG.data(c,:,d))', linspace(-1,1.999,3000), linspace(-1,1.999,3000), 1:100, 1000, 0.2);
                      %[tmp]=sub_mwt(squeeze(EEG.data(c,:,d))', linspace(-1,1.999,3000), linspace(-1,1.999,3000), 1:100, 1000, 2,0.3);
                      %specdata(c,:,:,d)=abs(tmp');
%                    end
                    tmp=abs(awt_freqlist(EEG.data(c,:,d),1000,1:100));
                    specdata(c,:,:,d)=abs(awt_freqlist(EEG.data(c,:,d),1000,1:100));
               end
           end
            %
        lowscore=index(1:round(length(index)/2));
        highscore=index((length(index)-round(length(index)/2))+1:end);
        lowindex=lowscore;highindex=highscore;
        if contains(preprocess,'hightrial')
            lowindex=highscore(1:round(length(highscore)/2));
            highindex=highscore(length(highscore)-round(length(highscore)/2)+1:end); 
        end
        specdata=permute(specdata,[2,3,1,4]);
        lfpt=linspace(-1,1.999,3000);
        %specdata=basecorrect(specdata,lfpt,-0.5,-0.1,'Subtract');
        specdata=permute(specdata,[3,1,2,4]);
        Highdata=squeeze(mean(specdata(:,:,:,highindex),4));
        Lowdata=squeeze(mean(specdata(:,:,:,lowindex),4));
end
function [Fvalue,Pvalue]=timepointanova(datamat)
for i=1:size(datamat,3)% for each channel
    parfor j=1:size(datamat,4) %for each time point
        tmpdata=squeeze(datamat(:,:,i,j,:));
        tmpdata=permute(tmpdata,[3,1,2]);
        anovaresult=simple_mixed_anova(tmpdata,[],{'location','intensity'});
        F_location(i,j)=anovaresult{3,4};
        F_intensity(i,j)=anovaresult{5,4};
        F_interaction(i,j)=anovaresult{7,4};
        P_location(i,j)=anovaresult{3,5};
        P_intensity(i,j)=anovaresult{5,5};
        P_interaction(i,j)=anovaresult{7,5};
    end
end
Fvalue{1}=F_location; Fvalue{2}=F_intensity; Fvalue{3}=F_interaction;
Pvalue{1}=P_location;Pvalue{2}=P_intensity;Pvalue{3}=P_interaction;
end
function data=evalvar(datamat,varname)
    data=eval(['datamat.',varname]);
end
function Savevariable(datamat,varname,var)
for i=1:length(varname)
    eval(['datamat.',varname{i},'=var{i};']);
end
end
function bool=Checkvariable(datafilename,varname)
    bool=false;
    if ~exist(datafilename,'file')
        bool=true;
    else
        tmpmat=matfile(datafilename);
        varnames=fieldnames(tmpmat);
        v=ismember(varname,varnames);
        if prod(v)==0
            bool=true;
        end
    end
end
%% auto remove ICs using IClabel 
% ...remove, 
function ICAremoval
subject=[2:74,76:109];
Location={'RH','LH'};
basepath='/mnt/Share/yuelp/Data_Fourmodalities/ICA';
savepath='/mnt/Share/yuelp/Data_Fourmodalities/ICAremovalnew'; 
mkdir(savepath);
parfor i=1:length(subject)
    EEG=ICAremove(strcat('S',num2str(subject(i)),'L_ICA.set'),basepath);
    EEG=pop_saveset(EEG,'filename',strcat('S',num2str(subject(i)),'L_ICA.set'),'filepath',savepath);
    EEG=ICAremove(strcat('S',num2str(subject(i)),'R_ICA.set'),basepath);
    EEG=pop_saveset(EEG,'filename',strcat('S',num2str(subject(i)),'R_ICA.set'),'filepath',savepath);
end
end
function EEGnew=ICAremove(filename,basepath)
    EEG = pop_loadset('filename',filename,'filepath',basepath);
    EEG = eeg_checkset( EEG );
    EEG = pop_iclabel(EEG, 'default');
    EEG = eeg_checkset( EEG );
    EEG = pop_icflag(EEG, [0 0.1;0.6 1;0.8 1;NaN NaN;NaN NaN;0.8 1;NaN NaN]);
    % muscle prob 0.6-1, eye prob 0.8-1 channel noise 0.8-1
    EEG = eeg_checkset( EEG );
    EEGnew=pop_subcomp(EEG,find(EEG.reject.gcompreject==1),0);
    EEGnew=eeg_checkset(EEGnew);
end
function [ERPHighdata,ERPLowdata,SpecHighdata,SpecLowdata,chanlocs]=getERPandSpecdata(filename,basepath,preprocess,trialaverage,specprocess,freqaverage,calspec)
       % get original EEG data from each subject- no filter
       % reassign trials to high and low level according to the pain
       % rating.
       % no baseline correct for TCA decomposition, baselinecorrect for
       % general ANOVA analysis.
       EEG_all = pop_loadset('filename',filename,'filepath',basepath);
       chanlocs=EEG_all.chanlocs;
       if contains(preprocess,'Fz')
       EEG_all = pop_reref(EEG_all, 'Fz');
        EEG_all = eeg_checkset( EEG_all);
       elseif contains(preprocess,'average')
           EEG_all=pop_reref(EEG_all,[]);
           EEG_all=eeg_checkset(EEG_all);
       end
       Event=struct2table(EEG_all.event);
       eventdescription=Event.type;  
       Event=struct2table(EEG_all.event);
       eventdescription=Event.type;
       modalityall={'L','S','V','A'};
       for mod=1:4
           modality=modalityall{mod};
       modalityindex=contains(eventdescription,modality);
       score=Event.rating;
       invalid=isnan(score(:,1));
       modalityindex(invalid)=false;
       EEG=pop_select(EEG_all,'trial',find(modalityindex==1));
       Event=struct2table(EEG.event);
       score=Event.rating;
       score=score(:,1);
        [~,index]=sort(score);
        eventdescription(isnan(score))=[];
        score(isnan(score))=[];
        EEG.data(:,:,isnan(score))=[];
        if mod==2
            EEG.data=arfit2interpolate(EEG.data,[950 1030],5);
        end
        lfpt=linspace(-1,1.999,3000);
        if calspec
        specdata=zeros(size(EEG.data,1),3000,100,size(EEG.data,3)); % calculate spectrogram
           for c=1:size(EEG.data,1)
               parfor d=1:size(EEG.data,3)
                   if contains(specprocess,'stft')
                      [~, tmp] = sub_stft(squeeze(EEG.data(c,:,d))', linspace(-1,1.999,3000), linspace(-1,1.999,3000), 1:100, 1000, 0.2);
                   elseif contains(specprocess,'mwt')
                      [tmp]=sub_mwt(squeeze(EEG.data(c,:,d))', linspace(-1,1.999,3000), linspace(-1,1.999,3000), 1:100, 1000, 2,0.3);
                      tmp=(tmp.*conj(tmp))/1000;
                   end
                      specdata(c,:,:,d)=abs(tmp');
               end
           end
        end
        EEGdata=permute(EEG.data,[2,1,3]);
        
        [~,index]=sort(score);
        lowscore=index(1:round(length(index)/2));
        highscore=index((length(index)-round(length(index)/2))+1:end);
        lowindex=lowscore;highindex=highscore;
%         HighScore=score(highindex);LowScore=score(lowindex);
        ERPHighdata(:,:,mod)=mean(EEGdata(:,:,highindex),3);
        ERPLowdata(:,:,mod)=mean(EEGdata(:,:,lowindex),3);
        if calspec
            specdata=permute(specdata,[2,3,1,4]);
            specdata=averagefreq(specdata);
            SpecHighdata(:,:,:,mod)=mean(specdata(:,:,:,highindex),4);
            SpecLowdata(:,:,:,mod)=mean(specdata(:,:,:,lowindex),4);
        end
        end
end
function specdatanew=averagefreq(specdata)
    % average frequnecy [1,4],[4,8],[8,13],[13,30],[30,50],[50,100]
    % for reduced data load in TCA decomposition
    freq={[1,4],[4,8],[8,13],[13,30],[30,50],[50,100]};
    frequency=linspace(1,100,size(specdata,2));
    for i=1:length(freq)
        specdatanew(:,i,:,:)=mean(specdata(:,frequency>=freq{i}(1)&frequency<=freq{i}(2),:,:),2);
    end
end
