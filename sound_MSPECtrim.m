function [mspec,f]=sound_MSPECtrim(y,durseg,numseg,fs,win,ovlp)
% 
% The function sound_MSPECtrim returns the mean spectra of the
% pressure corrected (uPa) sound data in y. y first is divided into 
% DURSEG (e.g., 5 sec) long non overlapping segments, and the NUMSEG 
% (e.g., 8) with the lowest broadband RMS are retained 
%  and used to calculate the mean spectra.  
% 
% EXAMPLE: 
% [mspec,f]=sound_MSPECtrim(y,5, 8,fs,2^14,0);
% 
% INPUT 
% y = pressure corrected time series uPa (e.g. 2 minute recording) 
% fs = sample rate Hz (default: 48000 Hz) 
%  
% Segment parameters: 
%   durseg   = duration of data segment to consider 
%              default = 5 seconds 
%   numseg   = number of segments(with ranked lowest rms)to retain 
%              default = 8.  If data are not long enough to keep 
%              the number requested, all the segments are retained 
%
% FFT parameters:  
%   win = number of data points in each FFT window (or NFFT). Should be
%        a power of 2. Default is 2^14 points  --> 2.93 Hz freq.
%        resolution at fs=48000 
%    ovlp = # points of overlap in FFT calcs; must be < win  
%        typically = 0 
% 
% OUTPUT 
%  mspec = mean power spectrum scaled to preserve mean square pressure   
%          Use 10*log10(mspec) to convert to dB. 
%          
%  f = frequency bins corresponding to mspec 
% 
% NOTES:
%  1. A Hanning window is applied to each win before FFT calculation 
%  2. The spectrum is scaled such that the sum of the spectral bins is 
%    equal to the mean squared amplitude (variance) in the timeseries  
%    (Eq 6.7 Glover et al. Modeling Methods in Marine Science).  
%    So one can sum the mspec over a range of frequencies to get an 
%    estimate of the sound pressure level in that band.  
%
%
%    For example to find SPL in the 100-20000 Hz band: 
%    a=find(f > 100 & f < 20000); 
%    In decibels: SPL = 10*log10(sum(mspec(a))); 


%% Del Bohnenstiehl - NCSU 
%  Modified June 2016  
%  Modified Feb 2021 
%  drbohnen@ncsu.edu 

% comparisons with pwelch  
%[mspec,f]=sound_MSPEC(y,48000,2^14,0);
% [Pxx,f] = pwelch(y,hanning(2^14),0,2^14,48000,'onesided','power');
% [PSD,f] = pwelch(y,hanning(2^14),0,2^14,48000,'onesided','psd');
% freq resolution = 48000/2^14  
% mode(mspec./PSD) = freq resolution (2.9297) 
% mode(Pxx./PSD) = 1.5 times the frequency resolution (4.3943) 
% mode(Pxx./mspec) = 1.5  (so a slight diffence in linear unit for pwelch).
% 


%% Set Defaults 

if nargin < 6
  warning('not enough arguemnts') 
   return 
end
if isempty(durseg); durseg=5;  
    disp('Using 5 second segments');   end 
if isempty(numseg); numseg=8;  d
    disp('Keeping 8 segments'); end 
if isempty(fs); fs=48000; disp('Using 48000 Fs');  end 
if isempty(win); win= 2^14; disp('Using NFFT length 2^14'); end 
if isempty(ovlp); ovlp=0; disp('no overlap');end 

%% Break the timeseries into nfft WIN size segments & store in columns 
y=y-mean(y);                     % demean the time series 
x=buffer(y,win,ovlp,'nodelay');  % buffer in a matrix  
if x(end)==0; x(:,end)=[]; end   % if last colum is zero padded delete. 
[r,Nwin]=size((x));              % find size of matrix 
msub=repmat(mean(x,1),r,1);      % calculate the mean of each column 
x=x-msub;                        % remove the mean of each column

%% Calc FFT on each segment  
wo=hanning(r);           % hanning window 
zo=x.*repmat(wo,1,Nwin); % apply hanning window  
nfft=2^nextpow2(win);    % next power of 2, although it should be already  
Y=fft(zo,nfft,1);        % two sided FFT opperating on each column 
po=2*abs(Y).^2/(nfft*sum(wo.^2)); % scale for PSD accounting for windowing 
po=po(1:ceil(nfft/2)+1,:); po(1)=0; % take first 1/2 & zero DC 
[prows,~] = size(po);               % # rows in po. 
m=0:1:prows-1; f=m*fs/nfft;         % define the frequency axis

%% Chunk spectral power matrix into DURSEG length blocks  
winpo=ceil(durseg/(win/fs));  % number of NFFT windows in durseg seconds 
[~,c]=size(po);               % number of windows (columns) in po matrix  
k=1:winpo:c-winpo;            % start position of each durseg segment  

%% Find NUMSEG with lowest broadband rms 
segrms=nan(1,length(k)); 
for q=1:length(k) 
 segrms(q)= sum(mean(po(:,k(q):k(q)+winpo),2));  % broadband rms   
end

%% Order the segment blocks smallest to largest 
[~,I]=sort(segrms);       % numseg with the smallest RMS first 

if numseg > length(I)
disp('* Not enough data to keep requested number of segments') 
disp('averaging all segments *') 
end

numseg=min([numseg,length(I)]);  % number of segments to keep 

%% Concatenate spectra for the quietest numseg 
po2=[];  
for z=1:numseg
po2=cat(2,po2,po(:,k(I(z)):k(I(z))+winpo));
end
mspec=mean(po2,2); % average spectra 


end  % function end 



