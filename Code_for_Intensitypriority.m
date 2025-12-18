% simplified edition for the formal analysis of Intensity priority paper.
clear all;
processtypeall={'nose','Fz','averageref'}; % three types of preprocess, '' means nose ref.
subject=[1:4 6:8 10:18 20:24 26:28 30:44 46:49 53:104];
g_subj=[3 5 6 9 12 15 16 18 19 23 24 26 28 29 31 32 36 39 40 45 46 49 50 51 52 53 54 61 64 65 67 68 72 73 75 77 78 79 81 82 83 87 89 90 91 93 95];
Location={'RH','LH'};
basepath='./'; % for EEG dataset1 95 subjects
lfpt=linspace(-1,1.999,3000);

for all=2:3
    processtype=processtypeall{all};% [],'_hightrial'.....
    datamat=matfile(['./trialaverage_',processtypeall{all},'.mat'],'Writable',true);
    EEG=pop_loadset('filename',strcat(num2str(subject(1)),'_LH.set'),'filepath',basepath);
    if all==2
        EEG=pop_reref(EEG,'Fz');
    end
    datamat.chanlocs=EEG.chanlocs;
    ERPsubject=nan(3000,length(EEG.chanlocs),4,95);Specsubject_mwt=nan(3000,6,length(EEG.chanlocs),4,95);Specsubject_stft=nan(3000,6,length(EEG.chanlocs),4,95);   
   parfor i=1:length(subject)
        [ERPLeftHigh,ERPLeftLow,SpecLeftHigh,SpecLeftLow,ScoreLeftHigh,ScoreLeftLow]=getERPandSpecdata(strcat(num2str(subject(i)),'_LH.set'),basepath,processtype,true,'mwt',true);
        [ERPRightHigh,ERPRightLow,SpecRightHigh,SpecRightLow,ScoreRightHigh,ScoreRightLow]=getERPandSpecdata(strcat(num2str(subject(i)),'_RH.set'),basepath,processtype,true,'mwt',true);
        ERPsubject(:,:,:,i)=cat(3,ERPLeftHigh,ERPRightHigh,ERPLeftLow,ERPRightLow);
        Specsubject_mwt(:,:,:,:,i)=cat(4,SpecLeftHigh,SpecRightHigh,SpecLeftLow,SpecRightLow);
        [~,~,SpecLeftHigh,SpecLeftLow]=getERPandSpecdata(strcat(num2str(subject(i)),'_LH.set'),basepath,processtype,true,'stft',true);
        [~,~,SpecRightHigh,SpecRightLow]=getERPandSpecdata(strcat(num2str(subject(i)),'_RH.set'),basepath,processtype,true,'stft',true);
        Specsubject_stft(:,:,:,:,i)=cat(4,SpecLeftHigh,SpecRightHigh,SpecLeftLow,SpecRightLow);
   end
datamat.ERPsubject=ERPsubject;
datamat.Specsubject_mwt=Specsubject_mwt;
datamat.Specsubject_stft=Specsubject_stft;
EEG=pop_loadset('filename',strcat(num2str(subject(1)),'_LH.set'),'filepath',basepath);
datamat.chanlocs=EEG.chanlocs;
end
%% get the tensormatrix (combined the condition and subject dimension)
for all=1:2
    processtype=processtypeall{all};% [],'_hightrial'.....
    datamat=matfile(['./trialaverage_',processtypeall{all},'.mat'],'Writable',true);
    spec=datamat.Specsubject_mwt;
    tensormatrix=[];
    for i=1:4
        tensormatrix=cat(4,tensormatrix,squeeze(spec(:,:,:,i,:)));
    end
    datamat.tensormatrix_mwt=tensormatrix;
    spec=datamat.Specsubject_stft;
    tensormatrix=[];
    for i=1:4
        tensormatrix=cat(4,tensormatrix,squeeze(spec(:,:,:,i,:)));
    end
    datamat.tensormatrix_stft=tensormatrix; % note that now the tiral order is LH, RH, LL, RL, different from previous data
end
%% Time anova for trialaveraged data
for all=1:2
datamat=matfile(['./trialaverage_',processtypeall{all},'.mat'],'Writable',true);
% condition order: LeftHigh,RightHigh,LeftLow,RightLow
Intensity={'High','Low'};
Side={'Left','Right'};
ERP=[];
ERPdata=datamat.ERPsubject;
ERP(1,1,:,:,:)=squeeze(ERPdata(:,:,1,:));
ERP(1,2,:,:,:)=squeeze(ERPdata(:,:,2,:));
ERP(2,1,:,:,:)=squeeze(ERPdata(:,:,3,:));
ERP(2,2,:,:,:)=squeeze(ERPdata(:,:,4,:));
[Fvalue,Pvalue]=timepointanova(ERP);
% for generating distribution, using resampled data for 20 times
F_resample=cell(20,3);
P_resample=cell(20,3);
parfor i=1:20
    [F_resample(i,:),P_resample(i,:)]=timepointanova(datasample(ERP,size(ERP,5),5));
end
datamat.ERP_Fvalue=Fvalue;
datamat.ERP_Pvalue=Pvalue;
datamat.ERP_Fvalue_resample=F_resample;
datamat.ERP_Pvalue_resample=P_resample;
end
%% get the significant latency of F value
for all=1:2
    datamat=matfile(fullfile('./',['trialaverage_',processtypeall{all},'.mat']));
    ERP_Fvalue=datamat.ERP_Fvalue;
    ERP_Fvalue_shuffle=datamat.ERP_Fvalue_resample;
    ERP_Pvalue=datamat.ERP_Pvalue;
    ERP_Pvalue_shuffle=datamat.ERP_Pvalue_resample;
    Pvaluecorrect=cellfun(@(x) pcorrect(x),ERP_Pvalue,'UniformOutput',0);
    Pvaluecorrect_shuffle=cellfun(@(x) pcorrect(x),ERP_Pvalue_shuffle,'UniformOutput',0);
    chanlocs=datamat.chanlocs;
    chanlocs=struct2table(chanlocs);
    chanlocs=chanlocs.labels;
    
    if all==1
       Fvalue=cellfun(@(x) x(:,ismember(chanlocs,'Cz')),ERP_Fvalue,'UniformOutput',0);
       Pvalue=cellfun(@(x) x(:,ismember(chanlocs,'Cz')),Pvaluecorrect,'UniformOutput',0);
       %Fvalue_shuffle=cellfun(@(x) x(:,ismember(chanlocs,'Cz')),ERP_Fvalue_shuffle,'UniformOutput',0);
       Pvalue_shuffle=cellfun(@(x) x(:,ismember(chanlocs,'Cz')),Pvaluecorrect_shuffle,'UniformOutput',0);
       sigtime_Cz=cellfun(@(x) getSignificant(x,lfpt),Pvalue,'UniformOutput',0);
       sigtime_Cz_shuffle=cellfun(@(x) getSignificant(x,lfpt),Pvalue_shuffle,'UniformOutput',0);
    elseif all==2
        Fvalue=cellfun(@(x) x(:,ismember(chanlocs,{'C3','C4'})),ERP_Fvalue,'UniformOutput',0);
        Pvalue=cellfun(@(x) x(:,ismember(chanlocs,{'C3','C4'})),Pvaluecorrect,'UniformOutput',0);
    end
    figure;
    subplot(1,2,1);plot(lfpt,Fvalue);
end

function p_correct=pcorrect(Pvalue)
    p_correct=reshape(mafdr(Pvalue(:),'BHFDR',true),size(Pvalue));
end
function [Fvalue,Pvalue]=timepointanova(datamat)
for i=1:size(datamat,3)% for each channel
    parfor j=1:size(datamat,4) %for each time point
        tmpdata=squeeze(datamat(:,:,i,j,:));
        tmpdata=permute(tmpdata,[3,1,2]);
        anovaresult=simple_mixed_anova(tmpdata,[],{'intensity','location'});
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
function [ERPHighdata,ERPLowdata,SpecHighdata,SpecLowdata,HighScore,LowScore,chanlocs]=getERPandSpecdata(filename,basepath,preprocess,trialaverage,specprocess,freqaverage)
       % get original EEG data from each subject- no filter
       % reassign trials to high and low level according to the pain
       % rating.
       % no baseline correct for TCA decomposition, baselinecorrect for
       % general ANOVA analysis.
       EEG = pop_loadset('filename',filename,'filepath',basepath);
       chanlocs=EEG.chanlocs;
       if contains(preprocess,'Fz')
        EEG = pop_reref( EEG, 'Fz');
        EEG = eeg_checkset( EEG );
       elseif contains(preprocess,'average')
           EEG=pop_reref(EEG,[]);
           EEG=eeg_checkset(EEG);
       end
       Event=struct2table(EEG.event);
       eventdescription=Event.type;  
       invalid=[];
       try
       invalid=cellfun(@(x) isempty(x),Event.rating,'UniformOutput',1); 
       score=cell2mat(Event.rating(~invalid));
       eventdescription=Event.type;
       eventdescription(invalid)=[];
       catch
          score=Event.rating; 
       end
        eventdescription(isnan(score))=[];
        score(isnan(score))=[];
        EEG.data(:,:,isnan(score))=[];
        lfpt=linspace(-1,1.999,3000);
       
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
%                     tmp=abs(awt_freqlist(EEG.data(c,:,d),1000,0:2:100));
%                    specdata(c,:,:,d)=abs(awt_freqlist(EEG.data(c,:,d),1000,0:2:100));
               end
           end
        EEGdata=permute(EEG.data,[2,1,3]);
        %EEGdata=basecorrect(EEGdata,lfpt,-0.5,0,'Subtract');
        specdata=permute(specdata,[2,3,1,4]);
        %specdata=basecorrect(specdata,lfpt,-0.5,-0.1,'Subtract');
        specdata=averagefreq(specdata);
        [~,index]=sort(score);
        lowscore=index(1:round(length(index)/2));
        highscore=index((length(index)-round(length(index)/2))+1:end);
        lowindex=lowscore;highindex=highscore;
        if contains(preprocess,'hightrial')
            lowindex=highscore(1:round(length(highscore)/2));
            highindex=highscore(length(highscore)-round(length(highscore)/2)+1:end); 
        end
        HighScore=score(highindex);LowScore=score(lowindex);
        if trialaverage
        ERPHighdata=mean(EEGdata(:,:,highindex),3);
        ERPLowdata=mean(EEGdata(:,:,lowindex),3);
        SpecHighdata=mean(specdata(:,:,:,highindex),4);
        SpecLowdata=mean(specdata(:,:,:,lowindex),4);
        else
            Highdata=EEGdata(:,:,highindex);
            Lowdata=EEGdata(:,:,lowindex);
            SpecHighdata=specdata(:,:,:,highindex);
            SpecLowdata=specdata(:,:,:,lowindex);
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
function sigtime=getSignificant(pvalue,t)
   %pvalue=mafdr(pvalue,'BHFDR',true);
    sigtime=t>0.1&t&pvalue'<0.05;
    c=1;siglength=0;sigt{1}=[];
    for i=1:length(sigtime)
        if sigtime(i)>0
            siglength(c)=siglength(c)+1;
            sigt{c}=cat(1,sigt{c},t(i));
        elseif sigtime(i)==0&&siglength(c)~=0
            c=c+1;
            siglength(c)=0;
            sigt{c}=[];
        end
    end       
        sigtime=min(cellfun(@(x) min(x),sigt(siglength>10),'UniformOutput',1));
        if isempty(sigtime)
           sigtime=nan;
        end
end

