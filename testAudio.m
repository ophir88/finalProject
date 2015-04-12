function testAudio(filename)

[mySignal,Fs] = audioread(filename);
signalTime = length(mySignal)/Fs;
L = length(mySignal)
windowSizeTime = 0.5; %seconds
overlapTime = 0.01;
overlapSize = overlapTime*Fs;
windowSize = windowSizeTime*Fs
A=[1,1];

for j=1:windowSize:L

NFFT = 2^nextpow2(L); % Next power of 2 from length of y
Y = fft(mySignal(j:j+windowSize),NFFT);
f = Fs/2*linspace(0,1,NFFT/2+1);
Y = Y.*conj(Y)/L;
% Plot single-sided amplitude spectrum.
subplot(4,4,round(j/windowSize+1))
plot(f,Y(1:NFFT/2+1)) 
title(sprintf( 'Single-Sided Amplitude Spectrum of y(t) %d', round(j/windowSize+1  )))
xlabel('Frequency (Hz)')
ylabel('|Y(f)|')
tranCoef = length(Y)/Fs;

minStart=max(Y)/35;
%minStart = 0.5;
inPeak = 0;
tempMax = 0;
currPeak = 0;
peaks = zeros(4,2);
peak_indx = 1;
 for i=round(27*tranCoef):round(4180*tranCoef)
 
    
    if(Y(i)>minStart)
        inPeak=1;
        
    else
        inPeak = 0;
    end
    
    if (inPeak==1)
        
        if(i>currPeak)
        k = i;
        while(Y(k)>minStart)
            if (Y(k)>tempMax)
            tempMax = Y(k);
            loc=k/tranCoef;
            end
            k=k+1;
        end
        peaks(peak_indx,1)=tempMax;
        peaks(peak_indx,2)=loc;
        peak_indx = peak_indx+1;
        if (peak_indx == 4)
            break;
        end
        currPeak=k;
        tempMax=0;
        end
        
        
    end
    
 
 end
 
 
 %look for 'side' peaks:
  for l=1:length(peaks)
     
     for m=1:length(peaks)
         if (m==l) 
             continue
         else
         if (peaks (l,2)>0 && abs(peaks(l,2)-peaks(m,2))<6 && peaks(l,1)-peaks(m,1)>peaks(l,1)/10)
             'found "side" peak'
             [peaks(m,1),peaks(m,2)]
             'and';
             [peaks(l,1),peaks(l,2)]
             peaks(m,2)=0;
         end
         end
         
         
     end
     
 end

 
 % look for octaves:
 for l=1:length(peaks)
     
     for m=1:length(peaks)
         if (m==l) 
             continue
         else
         if ((peaks (l,2)>0 && peaks(l,2)>(peaks(m,2)/2-2) && peaks(l,2)<(peaks(m,2)/2+2))||(peaks (l,2)>0 && peaks(l,2)>(peaks(m,2)/3-2) && peaks(l,2)<(peaks(m,2)/3+2)))
             'found duplicate'
             [peaks(m,1),peaks(m,2)]
             'and';
             [peaks(l,1),peaks(l,2)]
             peaks(m,2)=0;
         end
         end
         
         
     end
     
 end
 
 
A=vertcat(A,[9,9]);
 A=vertcat(A,[j/Fs,(j+windowSize)/Fs]);
 A=vertcat(A,[8,8]);
 A=vertcat(A,peaks);
    
    

end

A


% -----------------------------TESTING ZONE START:
% 
% for p=4:5
%     NFFT = 2^nextpow2(L); % Next power of 2 from length of y
% Y = fft(mySignal(windowSize*p:windowSize*(p+1)),NFFT);
% f = Fs/2*linspace(0,1,NFFT/2+1);
% % Y = .95.*Y./max(abs(Y));
% 
% Y = Y.*conj(Y)/L;
% % Plot single-sided amplitude spectrum.
% subplot(2,2,p-3)
% plot(f,Y(1:NFFT/2+1)) 
% title('Single-Sided Amplitude Spectrum of y(t)')
% xlabel('Frequency (Hz)')
% ylabel('|Y(f)|')
% tranCoef = length(Y)/Fs;
% minStart=max(Y)/50;
% %minStart = 0.5;
% inPeak = 0;
% tempMax = 0;
% currPeak = 0;
% peaks = zeros(4,2);
% peak_indx = 1;
%  for i=round(27*tranCoef):round(4180*tranCoef)
%  
%     
%     if(Y(i)>minStart)
%         inPeak=1;
%         
%     else
%         inPeak = 0;
%     end
%     
%     if (inPeak==1)
%         
%         if(i>currPeak)
%         k = i;
%         while(Y(k)>minStart)
%             if (Y(k)>tempMax)
%             tempMax = Y(k);
%             loc=k/tranCoef;
%             end
%             k=k+1;
%         end
%         peaks(peak_indx,1)=tempMax;
%         peaks(peak_indx,2)=loc;
%         peak_indx = peak_indx+1;
%         if (peak_indx == 4)
%             break;
%         end
%         currPeak=k;
%         tempMax=0;
%         end
%         
%         
%     end
%     
%  
%  end
%  'looking for dups in: '
%  peaks
% 
%  for l=1:length(peaks)
%      
%      for m=1:length(peaks)
%          if (m==l) 
%              continue
%          else
%          if (peaks (l,2)>0 && peaks(l,2)>(peaks(m,2)/2-2) && peaks(l,2)<(peaks(m,2)/2+2))
%              'found duplicate'
%              [peaks(m,1),peaks(m,2)]
%              'and'
%              [peaks(l,1),peaks(l,2)]
%              peaks(m,2)=0;
%          end
%          end
%          
%          
%      end
%      
%  end
%  peaks
% end
% 

 
    
    

%   ------------------------------------TEST END





's'
% [~,index] = max(y(1:m/2+1))
%     freq = 0:(Fs/m):Fs/2;
    % freq(index)
%     fprintf('Maximum occurs at %2.3f Hz\n',freq(index))




  


end