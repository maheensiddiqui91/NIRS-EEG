
function [v,stats] = GLM_general(EXPfiltdata_soc,m,chrom,ch,tt,rms_soc,new_hrf)
% This code has to be run for each fNIRS channel separately. 
% NirsSignal = fNIRS signal [samples x 1]
% X = design matrix [samples x (#regressors+1)]
% In this example, I consider the design matrix to be made of one regressor + costant term

% Generate random fNIRS signal at 5 Hz, 5 min long
NirsSignal=EXPfiltdata_soc(m).data(:,chrom,ch); %rand(1500,1);
fs=1;
% Generate design matrix
% Here as an example, I'm considering a protocol made of 30sec task blocks and 30 sec rest
stim=zeros(size(NirsSignal));

%tasktimes=30:60:260; % Onset of task blocks
%tasktimes=tt; %30:60:260; % Onset of task blocks
tasktimes2 = struct2cell(tt);
tasktimes = cell2mat(tasktimes2);
for i=1:length(tasktimes)
    stim( round(tasktimes(i)*fs): round((tasktimes(i)+12)*fs))=1;
end
%convolve with hrf
[hrf2] = spm_hrf(1/fs); %replace with EEG 
[hrf] = rms_soc; %rms2_soc_grand_avg(9:14,29,4);
Convolution1= conv(new_hrf,stim);
Convolution2 = conv(Convolution1,hrf);
%Convolution = conv(hrf,hrf2);
Convolution2= Convolution2(1:length(NirsSignal),:);
Convolution2=(Convolution2./max( Convolution2)); % 1 a.u.

%Design matrix
X=[Convolution2 ones(size(NirsSignal,1),1)];


nTRs=size(NirsSignal,1); % number of samples fNIRS signal
ndf=nTRs-rank(X); % Degrees of freedom
v=inv(X'*X)*X'*NirsSignal; % Compute all beta values for all regressors
stats=[];
stats.Var=(NirsSignal-X*v)'*(NirsSignal-X*v)/(nTRs-1-size(X,2)); % Compute variance
c=[1; 0]; % Define contrast, e.g. here I'm testing the significance of task vs rest
stats.t_stat=c'*v/sqrt(stats.Var*c'*inv(X'*X)*c); % Compute t-value from the betas
stats.p_value=1-tcdf(stats.t_stat,ndf); % Compute corresponding p-value

%figure; plot(Convolution2); hold on; plot(NirsSignal);
end

