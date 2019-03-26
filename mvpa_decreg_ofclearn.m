% mvpa_decreg_ofclearn.m
% 
% Wimmer & Büchel bioRxiv 2018 ROI mvpa code

% download and unzip decval_roifiles.zip
% download and unzip regfiles.zip


clear

% set start subject and end subject
subjstart =  2; % 2
subjend =   38; % 38
% set excluded subjects (movement, sleeping, chance performance, repeat task performance
exclude = [1 30 35 36];
exclude = [exclude 11 15 28]; % 11, 15, 28 motion

% set paths
pathmain = pwd;
pathreg = fullfile(pathmain, 'regfiles');
pathdec = fullfile(pathmain, 'decval_roifiles');

% set a few variables
nseq = 8; % 8 maze sequences
nrep = 4; % 4 repetitions of each maze sequence
% in case of very extreme values, censor
exclthresh = 5;
% do exclusion based on below-zero performance on actual state 1 and state 3 stimuli
doexclude = 1;
if doexclude
    fprintf(['\n' 'DOING base classifier performance exclusion' '\n'])
else
    fprintf(['\n' 'no base classifier performance exclusion' '\n'])
end


% primary ROIs
roistruct{1} = 'ofcc10mm'; % ofc 10mm
roistruct{2} = 'aalhipp'; % aal hippocampus
roistruct{3} = 'threecomb'; % occipital & temporal visual face, scene, object conjunction
roistruct{4} = 'hipp_ho25'; % harvard-oxford hippocampus 25% threshold

% 25 control PFC regions
roistruct{5} = 'c_rantofc';
roistruct{6} = 'c_lantofc';
roistruct{7} = 'c_gACC';
roistruct{8} = 'c_dfpc';
roistruct{9} = 'c_fpc';
roistruct{10} = 'c_rlatPFC';
roistruct{11} = 'c_llatPFC';
roistruct{12} = 'c_rlatOFC';
roistruct{13} = 'c_llatOFC';
roistruct{14} = 'c_ACC';
roistruct{15} = 'c_sACCx';
roistruct{16} = 'c_rmidPFC';
roistruct{17} = 'c_lmidPFC';
roistruct{18} = 'c_rinfPFC';
roistruct{19} = 'c_linfPFC';
roistruct{20} = 'c_ACCb';
roistruct{21} = 'c_ACCc';
roistruct{22} = 'c_rdlatPFC';
roistruct{23} = 'c_ldlatPFC';
roistruct{24} = 'c_rdlPFCa';
roistruct{25} = 'c_ldlPFCa';
roistruct{26} = 'c_rdlPFCb';
roistruct{27} = 'c_ldlPFCb';
roistruct{28} = 'c_rdlPFCc';
roistruct{29} = 'c_ldlPFCc';

% 4 control hippocampal subregions
roistruct{30} = 'rhippant20_ho25';
roistruct{31} = 'lhippant20_ho25';
roistruct{32} = 'rhipppost20_ho25';
roistruct{33} = 'lhipppost20_ho25';



% set roi here (or edit to loop across all)
roi = roistruct{1};


% load roi file
cd(pathdec)
load(['decval_' roi])
mvpadecval = decvalstruct;


% loop across category contrasts and subjects
for iRoi = 1:3
    
    if iRoi==1
        stimcontrast = 'fvs_'; stimfocus = 1; % face v scene
        decvalall = mvpadecval.fvs;
    elseif iRoi==2
        stimcontrast = 'svf_'; stimfocus = 2; % scene v face
        decvalall = mvpadecval.svf;
    elseif iRoi==3
        stimcontrast = 'ovf_'; stimfocus = 3; % obj v face
        decvalall = mvpadecval.ovf;
    end
    
    fprintf(['\n' stimcontrast '\n']);
    
    l = 1;
    
    for subj = subjstart:subjend
    if ismember(subj,exclude)
    else
        
        fprintf(['\n' 'subj = ' num2str(subj) '\n']);
        
        % set decvalues for this subject
        decval = decvalall{subj};
        
        % decision value: exclude extreme outliers >5 std from mean
        decval(decval>(nanmean(decval)+nanstd(decval)*exclthresh) | decval<(nanmean(decval)-nanstd(decval)*exclthresh)) = NaN;
        nexcludecheck = sum(isnan(decval))
        
        % invert decval variable for later
        negmult = -1;
        
        
        % load regressor information
        cd(pathreg)
        if subj<10; zeroa='0'; else; zeroa=''; end
        regind=load(['regs_learn_s' zeroa num2str(subj) '.mat']);
        cd(pathdec)
        
        % set regressors info
        stimlist = regind.stimlist;
        reward = regind.reward;
        stimstate1 = regind.stimstate1;
        stimstate3 = regind.stimstate3;
        stimnumber = regind.stimnumber;
        learnedorig = regind.learnedorig;
        
        reward(reward==0) = NaN;
        rewardorig = reward+(1-max(reward));
        reward0 = (rewardorig-0.5).*2;
        
        % get stim category directly from regressor file
        stim1face = stimstate1==1;
        stim1scene = stimstate1==2;
        stim1obj = stimstate1==3;
        
        stim3face = stimstate3==1;
        stim3scene = stimstate3==2;
        stim3obj = stimstate3==3;
        
        % pick out 5 different time periods of each trial from full list
        stim1 = NaN(length(stim1face),1);
        m = 1;
        n = 1;
        o = 1;
        p = 1;
        q = 1;
        for i = 1:length(stimnumber)
            if rem(i,5)==1
                decval1(m,1) = decval(i);
                m = m+1;
            end
            if rem(i,5)==2
                decval2(n,1) = decval(i);
                n = n+1;
            end
            if rem(i,5)==3
                decval3(o,1) = decval(i);
                o = o+1;
            end
            if rem(i,5)==4
                decval4(p,1) = decval(i);
                p = p+1;
            end
            if rem(i,5)==0
                decval5(q,1) = decval(i);
                q = q+1;
            end
        end
        
        
        % set columns stim matrix
        colslseq = 1;
        colslstim = 2;
        colslrep = 15;
        
        % set columns regression matrix
        colmeanrew = 1;
        colstim = 2;
        colseq = 3;
        colrep = 4;
        colrepcorr = 5;
        coltrial = 9;
        colrew = coltrial+1;
        coldecval = coltrial+2;
        colface1 = coldecval+1;
        colscene1 = coldecval+2;
        colobj1 = coldecval+3;
        colface3 = coldecval+4;
        colscene3 = coldecval+5;
        colobj3 = coldecval+6;
        colfacerew = coldecval+7;
        colscenerew = coldecval+8;
        colobjrew = coldecval+9;
        
        % set columns temp
        coltempseq = 2;
        colreptemprew = 5;
        
        excludetrials = find(isnan(stimlist(:,colslseq)));
        
        % add in sequence variable, repetition variable, reward, and decval variable.
        stim1(:,colmeanrew) = nanmean(reward0);
        stim1(:,colstim) = stimlist(:,colslstim);
        stim1(:,colseq) = stimlist(:,colslseq);
        stim1(:,colrep) = stimlist(:,colslrep);
        stim1(:,colrew) = reward0;
        stim1(:,coldecval) = decval1(:,1).*negmult;
        stim1(:,coltrial) = 1:40;
        
        
        repcorr = [];
        temp = [(1:size(stimlist,1))' stimlist(:,colslseq) stimlist(:,colslrep) learnedorig rewardorig];
        temp = sortrows(temp,coltempseq); % 2 = colseq
        
        for m = 1:nseq
            tempseq = find(temp(:,coltempseq)==m);
            tempseq = temp(tempseq,:);
            reptemp = zeros(nrep,1);
            for n = 2:nrep
                if tempseq(n-1,colreptemprew)==1
                    reptemp(n) = reptemp(n-1)+1;
                elseif reptemp(n-1)>0
                    reptemp(n) = reptemp(n-1)+1;
                else
                    reptemp(n) = reptemp(n-1);
                end
            end
            repcorr = [repcorr; reptemp];
        end
        
        
        % add in nans for seq9 and seq10
        repcorr(end+1:end+8) = NaN;
        temp(:,end+1) = repcorr; % repcorr = 6
        temp = sortrows(temp,1); % resort by trial
        repcorr = temp(:,6); repcorro = repcorr;
        repcorr = repcorr-nanmean(repcorr); 
        
        % add to stim1 and stim3
        stim1(:,colrepcorr) = repcorr;
        
        % basic stimulus contrasts
        stim1(:,colface1) = stim1face-nanmean(stim1face);
        stim1(:,colscene1) = stim1scene-nanmean(stim1scene);
        stim1(:,colobj1) = stim1obj-nanmean(stim1obj);
        
        stim1(:,colface3) = stim3face-nanmean(stim3face);
        stim1(:,colscene3) = stim3scene-nanmean(stim3scene);
        stim1(:,colobj3) = stim3obj-nanmean(stim3obj);
        
        % get forward effects X reward
        stim1facerew = stim1(:,colface1); stim1facerew(reward<0) = NaN;
        stim1facerew = stim1facerew-nanmean(stim1facerew); stim1facerew(isnan(stim1facerew)) = 0;
        
        stim1scenerew = stim1(:,colscene1); stim1scenerew(reward<0) = NaN;
        stim1scenerew = stim1scenerew-nanmean(stim1scenerew); stim1scenerew(isnan(stim1scenerew)) = 0;
        
        stim1objrew = stim1(:,colobj1); stim1objrew(reward<0) = NaN;
        stim1objrew = stim1objrew-nanmean(stim1objrew); stim1objrew(isnan(stim1objrew)) = 0;
        
        % remove mean from repetition variable
        stim1(:,colrep) = stim1(:,colrep)-nanmean(stim1(:,colrep));
        
        % create stim3
        stim3 = stim1; stim3(:,colstim) = 0;
        stim3(:,coldecval) = decval3.*negmult;
        % create stim2
        stim2 = stim1; stim2(:,colstim) = 0;
        stim2(:,coldecval) = decval2.*negmult;
        
        %%% run within-participant models to evaluate consistency / find outliers
        if stimfocus==1
            colfocus = colface1;
            colalt = colscene1;
            colf3 = colface3;
            corrvar1 = stim1face-stim1scene;
            corrvar3 = stim3face-stim3scene;
            stim1tgt = stim1(:,colface1);
            stim3tgt = stim3(:,colface3);
        elseif stimfocus==2
            colfocus = colscene1;
            colalt = colface1;
            colf3 = colscene3;
            corrvar1 = stim1scene-stim1face;
            corrvar3 = stim3scene-stim3face;
            stim1tgt = stim1(:,colscene1);
            stim3tgt = stim3(:,colscene3);
        elseif stimfocus==3
            colfocus = colobj1;
            colalt = colface1;
            colf3 = colobj3;
            corrvar1 = stim1obj-stim1face;
            corrvar3 = stim3obj-stim3face;
            stim1tgt = stim1(:,colobj1);
            stim3tgt = stim3(:,colobj3);
        end
        
        % concatenate state1 and state3
        corrvar1and3 = [corrvar1; corrvar3];
        stim1stim3 = [stim1(:,1:colobj3); stim3(:,1:colobj3)];
        
        % check if dec values relate to category of presented stimulus
        [rcritical, ~] = corr(stim1stim3(:,coldecval),corrvar1and3,'rows','complete');
        zcritical = fisherz(rcritical);
        [coef, ~, ~] = glmfit(corrvar1and3,stim1stim3(:,coldecval));
        coef = coef';
        
        decvalcorr1and3(l,1:5) = [subj NaN NaN zcritical coef(2)];
        
        if zcritical<0 && doexclude
            perfexclude = ones(length(stim1(:,coldecval)),1);
        elseif doexclude
            perfexclude = zeros(length(stim1(:,coldecval)),1);
        end
        
        % after the above checks...
        % remove seq9 and seq10 from decval
        stim1(excludetrials,coldecval) = NaN;
        stim3(excludetrials,coldecval) = NaN;
        
        
        
        
        %%% bin decision values by repcorr given stimulus of interest
        % deccorr1mean = state 1 (stimulus on screen); deccorr3mean = state 3
        for m = 1:4
            deccorr1mean(l,m) = nanmean(stim1(stim1(:,colfocus)>0 & repcorro==m-1,coldecval));
            
            deccorr3mean(l,m) = nanmean(stim1(stim1(:,colf3)>0 & repcorro==m-1,coldecval));
        end
        
        if perfexclude & doexclude
            deccorr1mean(l,:) = NaN;
            deccorr3mean(l,:) = NaN;
        end
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% within-subject regression %%%
        rep = stim1(:,colrep);
        repcorr = stim1(:,colrepcorr);
        
        if stimfocus>10
            negdecval = 1;
        else
            negdecval = -1;
        end
        
        
        [coefbasic,dev,stats] = glmfit([stim1tgt],decval1.*negdecval);
        coefbasic = coefbasic';
        
        [coefs1xa,dev,stats] = glmfit([stim1tgt,stim3tgt,repcorr,stim1tgt.*repcorr,stim3tgt.*repcorr],decval1.*negdecval);
        coefs1xa = coefs1xa';
        
        [coefs1ctrla,dev,stats] = glmfit([stim1tgt,repcorr,stim1tgt.*repcorr],decval1.*negdecval);
        coefs1ctrla = coefs1ctrla';
        
        [coefs1ctrlb,dev,stats] = glmfit([stim3tgt,repcorr,stim3tgt.*repcorr],decval1.*negdecval);
        coefs1ctrlb = coefs1ctrlb';
        
        [coefs2xa,dev,stats] = glmfit([stim1tgt,stim3tgt,repcorr,stim1tgt.*repcorr,stim3tgt.*repcorr],decval2.*negdecval);
        coefs2xa = coefs2xa';
        
        [coefs3xa,dev,stats] = glmfit([stim1tgt,stim3tgt,repcorr,stim1tgt.*repcorr,stim3tgt.*repcorr],decval3.*negdecval);
        coefs3xa = coefs3xa';
        
        
        
        % exclude based on below-zero correlation between decision values and state 1 and state 3 category
        if perfexclude & doexclude
            coefbasic(:) = NaN;
            coefs1xa(:) = NaN;
            coefs1ctrla(:) = NaN;
            coefs1ctrlb(:) = NaN;
            coefs2xa(:) = NaN;
            coefs3xa(:) = NaN;
            zctrl = NaN;
        end
        
        decvalcoefbasic(l,1) = subj;
        decvalcoefbasic(l,2) = coefbasic(2);
        
        
        decvalcoefs1xa(l,1) = subj;
        decvalcoefs1xa(l,2:6) = coefs1xa(2:6);
        
        decvalcoefs1ctrla(l,1) = subj;
        decvalcoefs1ctrla(l,2:4) = coefs1ctrla(2:4);
        decvalcoefs1ctrlb(l,1) = subj;
        decvalcoefs1ctrlb(l,2:4) = coefs1ctrlb(2:4);
        
        decvalcoefs2xa(l,1) = subj;
        decvalcoefs2xa(l,2:6) = coefs2xa(2:6);
        decvalcoefs3xa(l,1) = subj;
        decvalcoefs3xa(l,2:6) = coefs3xa(2:6);
        
        
        clear coef*
        
        % concatenate
        subjnum = zeros(length(stim1),1); subjnum(:) = subj;
        stim1all = [subjnum perfexclude stim1];
        stim3all = [subjnum perfexclude stim3];
        stim2all = [subjnum perfexclude stim2];
        
        if exist('decval1all')
            decval1all = [decval1all; stim1all];
        else
            decval1all = [stim1all];
        end
        if exist('decval3all')
            decval3all = [decval3all; stim3all];
        else
            decval3all = [stim3all];
        end
        
        if exist('decval2all')
            decval2all = [decval2all; stim2all];
        else
            decval2all = [stim2all];
        end
        
        
        l = l+1;
    end
    end
    
    fprintf(['\n' stimcontrast '\n']);
    
    
    if strcmp(stimcontrast,'fvs_')
        contrast = 'decvalall_facevscene_';
    elseif strcmp(stimcontrast,'svf_')
        contrast = 'decvalall_scenevface_';
    elseif strcmp(stimcontrast,'ovf_')
        contrast = 'decvalall_objvface_';
    end
    
    cd(pathmain)
    cmd = ([' save ' contrast roi ' decvalcorr1and3 decval1all decval3all decval2all deccorr1mean deccorr3mean decvalcoef*']); eval(cmd);
    
    clear decval*
    
    fprintf(['\n' 'ROI: ' roi '\n'])
    if doexclude
        fprintf(['\n' 'DOING base classifier performance exclusion' '\n'])
    else
        fprintf(['\n' 'no base classifier performance exclusion' '\n'])
    end
    
end


% when finished, load and display some results



% extract values from three contrasts
cd(pathmain)
fvs = load(['decvalall_facevscene_' roi]);
svf = load(['decvalall_scenevface_' roi]);
ovf = load(['decvalall_objvface_' roi]);

coefs1xa(1,:,:) = fvs.decvalcoefs1xa(:,2:6);
coefs1xa(2,:,:) = svf.decvalcoefs1xa(:,2:6);
coefs1xa(3,:,:) = ovf.decvalcoefs1xa(:,2:6);

% average across three contrasts
for iS = 1:size(coefs1xa,2)
    decsubj = squeeze(coefs1xa(:,iS,:));
    coefs1xamean(iS,:) = nanmean(decsubj);
end

% coefs1xamean = key result matrix
nanmean(coefs1xamean)

% show ttest results for state 1 (even though this is biased by exclusion)
[h, p, ci, stats] = ttest(coefs1xamean(:,1));
stats.p = p;
stats

% show ttest results for future state by repetition correct interaction
[h, p, ci, stats] = ttest(coefs1xamean(:,5));
stats.p = p;
stats

% bar plot
bar(nanmean(coefs1xamean))



