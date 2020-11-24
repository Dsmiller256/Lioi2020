clear all; clc;
subjectlist={'sub-xp101','sub-xp102','sub-xp103','sub-xp104','sub-xp105','sub-xp106','sub-xp107','sub-xp108','sub-xp109','sub-xp110'};

fulldataset=[];
for sub=1:length(subjectlist)
    subject=subjectlist{sub};
    load(['D:' filesep 'Lioi2020' filesep 'Lioi2020' filesep 'rawdata' filesep 'derivatives' filesep subject filesep 'NF_eeg' filesep 'd_' subject '_task-eegNF_NFeeg_scores.mat'])
    filtered_nf=NF_eeg.nf(1:4:end);
    fulldataset=[fulldataset;filtered_nf];
    
    
    fig=figure;
    plot(filtered_nf);
    xlabel('Time(s)')
    ylabel('NF Score')
    ylim([0 1]);
    xlim([0 400]);

    cd ('D:\Lioi2020\Lioi2020\analysis\beh_analysis')
    savefig(fig,[subject '_NFscore.fig']);
    
    close all;
end

fig=figure;
plot(mean(fulldataset));
xlabel('Time(s))');
ylabel('NF Scoore');
ylim([0 0.6]);
xlim([0 400]);

savefig(fig,['All_NFscore.fig']);


