function main(filename)


% ---------------- variables: ---------------
% read file:
[mySignal,Fs] = audioread(filename);
% file length in seconds:
signalTime = length(mySignal)/Fs;
% file length vector
L = length(mySignal)

% window size, will eventualy be set by the midi file tempo variable:
%window size in seconds:
windowSizeTime = 0.5; %seconds
%window size vector:
windowSize = windowSizeTime*Fs

%   overlap: currently not used:
% overlapTime = 0.01;
% overlapSize = overlapTime*Fs;

% threshold
threshold = .1;



% results: [between time, to time, note, amp]
results = [0,0,0,0]
results1 = zeros(floor(L/windowSize)+1,4);

% keep temp peaks during 'for' loop: [between time, and time,freq, amp]
peaks = zeros(6,4);
peaksIdx = 1;
resultIdx = 1;
newFreq = 0;
% ------------- start of function: ---------------------

% break the input into 'window size' elements, and analyze separatly
%   iterate over input in 'window size' jumps
for k=1 : windowSize : L
    
    % get segment (check if not last one)
    if (k+windowSize>L)
        spectrum = fft(mySignal(k:L),8*L);
    else
        spectrum = fft(mySignal(k:k+windowSize),8*L);
    end
    spectrum = abs(spectrum(1:ceil(length(spectrum)/2)));
    
    % max(spectrum);
    filter = gaussfilt(128);
    spectrum = spectrum .^ 2;
                spectrum2 = conv(spectrum,filter);
    % regular fft:
    spectrum = spectrum/max(spectrum);
    % filtered fft:
                spectrum2 = spectrum2/max(spectrum2);
    % convert between vector location to real frequency:
    tranCoef = Fs/(length(spectrum)*2);
    %plot signals:
        figure(round(k/windowSize)+1)
        plot(linspace(0,Fs/2,length(spectrum)),spectrum)
        axis([0 1500 0 1])
        figure(round(k/windowSize)+11)
        plot(linspace(0,Fs/2,length(spectrum2)),spectrum2)
        axis([0 1500 0 1])
    
    % ------------------- find peaks: ----------------
    
    % look only at relevant frequencies (piano notes)
    
    for i=round(27*length(spectrum)/Fs):round(4180*length(spectrum)/Fs)
        
        % max point:
        if(spectrum(i)>threshold & spectrum(i)>spectrum(i-1)& spectrum(i)>spectrum(i+1))
           'segment peaks'
        peaks 
            
            % make sure peaks are not too close
           
            if(abs(i*tranCoef -peaks(peaksIdx,3)) > 10)
                
                % close peaks, look for stronger:
                % case new one is stronger, overwrite
                % older
%                 if(spectrum(i) > peaks(peaksIdx,4))
                    %new peak is close and much stronger, reaplace it with the
                    %former one:
                    newFreq = i*tranCoef;
                    
                    %check that it is not an harmonic of a already added note:
                    octave = 0;
                    for j=1:peaksIdx
                        
                        % if the new frequency is a multiplication of an
                        % existing one
                       
                        if(peaks(j,3)>0 &&(mod(newFreq,peaks(j,3)) < .03*newFreq || ...
                                peaks(j,3) - mod(newFreq,peaks(j,3))  < .03*newFreq))
                            octave = 1;
                            'found mult'
                                newFreq
                                peaks(j,3)
                            % If the magnitude of the peak is strong enough, declare that it is a different note instead of a loud harmonic
                            if(spectrum(i) > .3*peaks(j,4))
                                'this amp'
                                spectrum(i)
                                'bigger than'
                                .3*peaks(j,4)
                                'but big enough'
                                
                                
                                % add peak:
                                %                                                                 'before add'
                                %                                                                 peaks
                                %                                                                 'after add'
                                peaks(peaksIdx,1)=round(k/windowSize)*windowSizeTime;
                                peaks(peaksIdx,2)=round((k/windowSize)+1)*windowSizeTime;
                                peaks(peaksIdx,3)=i*tranCoef;
                                peaks(peaksIdx,4)=spectrum(i);
%                                 peaks(peaksIdx)=[round(k/windowSize)*windowSizeTime,round((k/windowSize)+1)*windowSizeTime,i*tranCoef, spectrum(i)];
%                                 peaks = vertcat(peaks,[round(k/windowSize)*windowSizeTime,round((k/windowSize)+1)*windowSizeTime,i*tranCoef, spectrum(i)]);
                                peaksIdx = peaksIdx +1;
                                'added'
                                break;
                            end
                        
                        else
                        end
                            
                        
                    end
                    if( octave == 0)
                             
                        
                         peaks(peaksIdx,1)=round(k/windowSize)*windowSizeTime;
                                peaks(peaksIdx,2)=round((k/windowSize)+1)*windowSizeTime;
                                peaks(peaksIdx,3)=i*tranCoef;
                                peaks(peaksIdx,4)=spectrum(i);
%                         peaks(peaksIdx)=[round(k/windowSize)*windowSizeTime,round((k/windowSize)+1)*windowSizeTime,i*tranCoef, spectrum(i)];

%                             peaks = vertcat(peaks,[i*tranCoef, spectrum(i)]);
                              peaksIdx = peaksIdx +1;
                    end
                else continue;
%                 end
            end
            
        end
    end
    % ------------------- end find peaks -------------
    
    
%     (horzcat(round(k/windowSize)*windowSizeTime,round((k+1)/windowSize)*windowSizeTime,peaks))
    results(resultIdx)=peaks;
    resultIdx = resultIdx+1;
    peaks = [0,0,0,0];
    peaksIdx=1;
end

[rows,columns] = size(results);
    if(columns >= 2)
        results = results(:,2:columns);
    end

results


end

function filter=gaussfilt(N)

% calculate alpha so the filter fills all N points
alpha=N;
first=-(1-N/2)*exp(-(1-N/2)^2/alpha);
count=0;
while first<.1*(-(1530/4000*N-N/2)*exp(-(1530/4000*N-N/2)^2/alpha))
    count=count+1;
    alpha=N*500*count;
    first=-(1-N/2)*exp(-(1-N/2)^2/alpha);
end


for n=1:N
    filter(n)=-(n-N/2)*exp(-(n-N/2)^2/alpha);   % d/dt of a gaussian
end
filter=filter/sum(abs(filter));     % normalization
end
