
%%
% FINAL parameters based on the GLM fitting to NEST data in script

%                                                     defaults
%                                                    (seconds)
%   p(1) - delay of response (relative to onset)         6
%   p(2) - delay of undershoot (relative to onset)      16
%   p(3) - dispersion of response                        1
%   p(4) - dispersion of undershoot                      1
%   p(5) - ratio of response to undershoot               6
%   p(6) - onset (seconds)                               0
%   p(7) - length of kernel (seconds)                   32
%


p_HbO2(1) = 12;
p_HbO2(2) = 17;
p_HbO2(3) = 1;
p_HbO2(4) = 1;
p_HbO2(5) = 2; 
p_HbO2(6) = 1;
p_HbO2(7) = 30;

p_HHb(1) = 13;
p_HHb(2) = 20;
p_HHb(3) = 1;
p_HHb(4) = 1;
p_HHb(5) = 2;
p_HHb(6) = 1;
p_HHb(7) = 30;

p_CCO(1) = 12;
p_CCO(2) = 21; 
p_CCO(3) = 1;
p_CCO(4) = 1;
p_CCO(5) = 2;
p_CCO(6) = 1;
p_CCO(7) = 30;

hrf_1=spm_hrf(1,p_HbO2);
hrf_2=spm_hrf(1,p_HHb);
hrf_3=spm_hrf(1,p_CCO);

new_hrf = [hrf_1 hrf_2 hrf_3];

figure; plot(hrf_1,'r'); hold on; plot(-hrf_2,'b'); plot(hrf_3,'g');
%%

%The following will load the correct data which is obtained from "spectral_NEST_FFT_FINAL_291121" where missing values are filled
%in and then baseline correction is performed. The correct segments are
%obtained for both social and non-social.
load('/Users/maheensiddiqui/Documents/GLM_general/EEG_RMS_NEW_291121_soc_nonsoc.mat');
load('/Users/maheensiddiqui/Documents/GLM_general/NEST_files/GLM_NEST_datastruct_soc.mat'); 


EEG_Ch_list =  [1:32];
EEG_list_names = {'P7';'P3';'Cz';'Pz';'P4';'TP8';'O1';...
    'O2';'T8';'F8';'C4';'F4';'Fp2';'Fz';'C3';'F3';'Fp1';'T7';...
    'F7';'Oz';'PO4';'FC6';'C2';'P10';'CP6';'PO8';'PO7';'CP5';'TP7';'FC5';'P9';'PO3'};

for freq = 1:5
for m=1:15
    
    tt = DMAT_soc(m);
    
for chrom=1:3
    for ch=1:19
        for j=1:32
            if isnan(rms2_soc_avg(8:16,j,m,freq))
                beta(:,chrom,ch,j,m,freq) = NaN;
                stat(m).t_stat(chrom,ch,j,freq) = NaN;
                stat(m).p_value(chrom,ch,j,freq) = NaN;
            else 
       [v(:,chrom,ch,j,freq),stats]= GLM_general(EXPfiltdata_soc,m,chrom,ch,tt,rms2_soc_avg(8:16,j,m,freq),new_hrf(:,chrom));
       beta(:,chrom,ch,j,m,freq) = v(1,chrom,ch,j,freq);
       stat(m).t_stat(chrom,ch,j,freq) = stats.t_stat;
       stat(m).p_value(chrom,ch,j,freq) = stats.p_value;
            end
        end
    end
end
end
end

for freq=1:5
for chrom=1:3
    for ch=1:19
        for ek=1:32
            [h(chrom,ch,ek,freq),p(chrom,ch,ek,freq),ci,stats]...
            = ttest(beta(:,chrom,ch,ek,:,freq),0);
        
        group_stats_soc(chrom,ch,ek,freq) = stats.tstat;
        
        end
    end
end
end

%p = squeeze(p);
for freq=1:5
%for chrom=1:3
    for ch=1:19
        for ek=1:32
            FDR(:,ch,ek,freq) = mafdr(p(:,ch,ek,freq),'BHFDR',true);
        end
    end
%end
end


FDR2 = FDR;
for freq=1:5
for chrom=1:3 %window
      for ch=1:19 
          for ek=1:32 %NIRS chans
              
            %[ind1, ind2] = find(FDR(chrom,ch,ek)>0.05);
            %FDR2(ind1,chrom,ch,ek) = NaN;
            if FDR(chrom,ch,ek,freq)>0.05
                FDR2(chrom,ch,ek,freq) = NaN;
            else
                FDR2(chrom,ch,ek,freq) = FDR(chrom,ch,ek,freq);
            end
          end
      end
end
end



p2 = p;
for freq=1:5
for chrom=1:3 %window
      for ch=1:19 
          for ek=1:32 %NIRS chans
              
            %[ind1, ind2] = find(FDR(chrom,ch,ek)>0.05);
            %FDR2(ind1,chrom,ch,ek) = NaN;
            if p(chrom,ch,ek,freq)>0.05
                p2(chrom,ch,ek,freq) = NaN;
            else
                p2(chrom,ch,ek,freq) = p(chrom,ch,ek,freq);
            end
          end
      end
end
end



HbO2 = squeeze(FDR2(1,:,:,:));
CCO = squeeze(FDR2(3,:,:,:));
HHb = squeeze(FDR2(2,:,:,:));
HbO2_p = squeeze(p2(1,:,:,:));
CCO_p = squeeze(p2(3,:,:,:));
HHb_p = squeeze(p2(2,:,:,:));
% Before FDR correction
CCO_tstat = zeros(19,32,4);
HbO2_tstat = zeros(19,32,4);
HHb_tstat = zeros(19,32,4);
group_stats_soc=squeeze(group_stats_soc);

for freq=1:5
for i=1:19
    for j=1:32
    if isnan(CCO_p(i,j,freq))
        CCO_tstat(i,j,freq) = NaN;
        
    else 
        CCO_tstat(i,j,freq) = group_stats_soc(3,i,j,freq); 
    end
    end
end
end

for freq=1:5
for i=1:19
    for j=1:32
    if isnan(HbO2_p(i,j,freq))
        HbO2_tstat(i,j,freq) = NaN;
        
    else 
        HbO2_tstat(i,j,freq) = group_stats_soc(1,i,j,freq); 
    end
    end
end
end

for freq=1:5
for i=1:19
    for j=1:32
    if isnan(HHb_p(i,j,freq))
        HHb_tstat(i,j,freq) = NaN;
        
    else 
        HHb_tstat(i,j,freq) = group_stats_soc(2,i,j,freq); 
    end
    end
end
end
% For FDR p-values

CCO_tstat_FDR = zeros(19,32,4);
HbO2_tstat_FDR = zeros(19,32,4);
HHb_tstat_FDR = zeros(19,32,4);
group_stats_soc=squeeze(group_stats_soc);

for freq=1:5
for i=1:19
    for j=1:32
    if isnan(CCO(i,j,freq))
        CCO_tstat_FDR(i,j,freq) = NaN;
        
    else 
        CCO_tstat_FDR(i,j,freq) = group_stats_soc(3,i,j,freq); 
    end
    end
end
end

for freq=1:5
for i=1:19
    for j=1:32
    if isnan(HbO2(i,j,freq))
        HbO2_tstat_FDR(i,j,freq) = NaN;
        
    else 
        HbO2_tstat_FDR(i,j,freq) = group_stats_soc(1,i,j,freq); 
    end
    end
end
end

for freq=1:5
for i=1:19
    for j=1:32
    if isnan(HHb(i,j,freq))
        HHb_tstat_FDR(i,j,freq) = NaN;
        
    else 
        HHb_tstat_FDR(i,j,freq) = group_stats_soc(2,i,j,freq); 
    end
    end
end
end
%

tstats_allchrom(:,:,1,:) = HbO2_tstat;
tstats_allchrom(:,:,2,:) = HHb_tstat;
tstats_allchrom(:,:,3,:) = CCO_tstat;
tstats_allchrom_FDR(:,:,1,:) = HbO2_tstat_FDR;
tstats_allchrom_FDR(:,:,2,:) = HHb_tstat_FDR;
tstats_allchrom_FDR(:,:,3,:) = CCO_tstat_FDR;

tstats_allchrom(6:8,:,:,:) = NaN;
tstats_allchrom(10,:,:,:) = NaN;
tstats_allchrom(19,:,:,:) = NaN;
tstats_allchrom_FDR(6:8,:,:,:) = NaN;
tstats_allchrom_FDR(10,:,:,:) = NaN;
tstats_allchrom_FDR(19,:,:,:) = NaN;
tstats_allchrom(:,13,:,:) = NaN;
tstats_allchrom(:,17,:,:) = NaN;
tstats_allchrom_FDR(:,13,:,:) = NaN;
tstats_allchrom_FDR(:,17,:,:) = NaN;
%
%%
% load('tstats_allchrom_soc.mat')
% load('tstats_allchrom_nonsoc.mat')
%%
%tstats_allchrom_FDR_soc = tstats_allchrom_FDR;
nonsig_chans{1} = [1 2 3 5 6 7 8 9 10 15 16 17 18];%[4 12 13 14];
nonsig_chans{2} = [1:10 13 15:17]; %[11 12 14 18];
nonsig_chans{3} = [1:3 5:10 15:17]; %[4 11 12 13 14 18];

tstats_allchrom_FDR_soc(nonsig_chans{1},:,1,:) = NaN;
tstats_allchrom_FDR_soc(nonsig_chans{2},:,2,:) = NaN;
tstats_allchrom_FDR_soc(nonsig_chans{3},:,3,:) = NaN;

tstats_allchrom_FDR_soc(:,[10 11 12 13 16 17 19 22],:,:) = NaN;
tstats_allchrom_FDR_soc(:,[10 11 12 13 16 17 19 22],:,:) = NaN;

EEG_Ch_list =  [1:32];
EEG_list_names = {'P7';'P3';'Cz';'Pz';'P4';'TP8';'O1';...
    'O2';'T8';'F8';'C4';'F4';'Fp2';'Fz';'C3';'F3';'Fp1';'T7';...
    'F7';'Oz';'PO4';'FC6';'C2';'P10';'CP6';'PO8';'PO7';'CP5';'TP7';'FC5';'P9';'PO3'};
chroms = {'HbO_2'; 'HHb'; 'oxCCO';};
freqs = {'Theta'; 'Alpha'; 'Beta'; 'Gamma'};
for freq=1:5
    figure
for i=1:3
    
subplot(2,2,i); 

c=imagesc(tstats_allchrom_FDR_soc(:,:,i,freq)); 
hold all;
set(c,'AlphaData',~isnan(tstats_allchrom_FDR_soc(:,:,i,freq)))
caxis([-4 4]); 
xlim([0.5 32.5]);
%colormap(redblue);

xticklabels = EEG_list_names;
xticks = 1:1:32;
yticks = 1:1:19;
set(gca, 'XTick', xticks,'XTickLabel', xticklabels, 'YTick', yticks,'FontSize',8)
title(chroms{i},'FontSize',15);
 
end
end
