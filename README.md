# cc_calculations

# Matlab Code

% reduction of the noise level from raw signal
% set the amplitudes within certain range (etc. 10uV--15uV) to 0

%%%%% import data file
% data(): .txt (ASCII),which is converted from .mcd through MC_DataTool,
% unit of data is: ms
% convertion of mcd file must follow the adchID order, click the button 21,
% 31, 41,51,...58,68,78 to put them into the 'Selection' window and click 'Save' 
% note: reference electrode (#15) does not need to be slected.

str1 = '***** sampling frequency of data: ';
str2 = 'Hz';
fre_power = [];


new_data = data;   % the unit for the first column of data() is : ms

L = length(data(:,1));
Fs = 1000/(data(2,1)-data(1,1));   % unit: Hz
channels_no = size(data,2)-1;
NFFT = 2^nextpow2(L); % Next power of 2 from length of y

fprintf('%s %f %s', str1, Fs, str2 );


for  id=1:channels_no       % run all the channels with signals
    
%----  reduction of the noise of the signal
    
%     for i = 1:L
%        if abs(new_data(i,id+1))<=10    % downsampling: 5khz and 2khz: <=15;  1khz: <=12; 200hz: <=10
%            new_data(i,id+1) = 0;
%        end
%     end


%----- calculation of the frequency change with time
     h=window('hamming',128);
     [B,f,t]=spectrogram(new_data(:,id+1),h,60,256,Fs);
    
    % figure
% refer to package bfilter2
% Bilateral filtering method was used for edge-preserving smoothing.
% B = bfilter2(A,W,SIGMA)   (value A must be closed interval [0,1])
     w     = 5;       % bilateral filter half-width
     sigma = [3 0.1]; % bilateral filter standard deviations
     B0=abs(B)./50000;  % divided by a large value( 10000) to make sure B0 within [0,1]
     C = bfilter2(B0,w,sigma);
%   subplot(2,1,2);
%   pcolor(t,f,C);
%   shading interp;
%   title('Frequency Time Distribution ')
%   xlabel('Time (S)')
%   ylabel('Frequency (Hz)')
%   ylim([0 10])
  
    
% -----  calculation of the power spectrum
     
  Y = fft(new_data(:,id+1),NFFT)/L;
  f1 = Fs/2*linspace(0,1,NFFT/2);
  Am = abs(Y(1:NFFT/2)); % Amplitude,single-sided (NFFT/2)
  Po = 2*Am.*Am;            % Power
  fre_power = cat(2,fre_power,Po);
  
% % Plot power spectrum.
%    subplot(3,1,3);
%   plot(f1,Po)   % plot either amplitude or power
%   title('Power Spectrum ')
%   xlabel('Frequency (Hz)')
%   ylabel('Power') 
%   xlim([0 10])
 
  % Plot signals.
 % subplot(2,1,1);
%   plot( new_data(:,1)./1000, new_data(:,id+1));
%   plot( data(:,1)./1000, data(:,id+1));
%   title( ' Signals ') % with downsampling of 1 KHz ') 
%   xlabel('Time (s)')
%   ylabel(' voltage (\muV)')
  
  
end

  fre_power = cat(2,f1',fre_power);
  fre_0_to_10hz = fre_power(1:657,:);
    
  [cc,prob] = corrcoef(new_data(:,2:end)) ; 
  top_cc = triu(cc,1) ;    % Kth diagonal: start from k==1 of upper triangular part of matrix (exclude k==0), the main diagonal is k==0
  [top_cc_row,top_cc_col,top_cc_value] = find(top_cc);
  nonzero_top_cc = cat( 2, top_cc_row, cat(2,top_cc_col,top_cc_value));
  
  figure
 %pcolor(cc)
 imagesc(cc)
 colormap(jet)
 title('correlation coefficients DIV14 33220 control') %(50\muM Nicotine) (Baseline1) ')
  xlabel('Electrode No.')
  ylabel('Electrode No.')
  
% saves the data and figure--should change name to desired path every time.  
  savefile = 'C:\Users\James\Desktop\Ziyue\Recordings 3-2018\Control\DIV14_33220_control\DS_DIV14_33220_control.mat'; 
  save(savefile)
  savefigure = 'C:\Users\James\Desktop\Ziyue\Recordings 3-2018\Control\DIV14_33220_control\DS_DIV14_33220_control.fig'; 
  savefig(savefigure)
