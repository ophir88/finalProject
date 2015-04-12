% Takes the times of note hits, along with the notes that are played in
% each window and the sampling frequency, and returns a matrix
% representation of the sheet music along with the tempo the sheet music
% should be played at.
%
% Inputs:
%   times - the times that notes were hit at
%   notes - a matrix of the notes that were hit in each window
%   fs    - the sampling frequency
%
% Outputs:
%   music - a matrix representation of the sheet music. This contains the
%           notes played and the lengths of the notes.
%   tempo - the tempo that the returned music should be played at.
function [music,tempo] = noteDuration(times,notes,fs)
durs = zeros(1,length(times)-1);
for n = 1:length(durs)
    durs(n)=times(n+1)-times(n);
end
durs;
dist=zeros(1,ceil(max(durs)/1000)*1000);
for n = 1:length(durs)
    dist(durs(n))=dist(durs(n))+1;
end

N = 60;
histdist=histc(durs,0:ceil(length(dist)/N):length(dist));

bar([0:ceil(length(dist)/N):length(dist)]/fs,histdist);
title('Histogram of Note Durations','fontsize',28);
xlabel('Duration (s)','fontsize',28);
ylabel('Number of occurences','fontsize',28);

[val,loc]=max(histdist);
loc
val
loc = (loc-.5)*length(dist)/N
quarter=avex2(dist(ceil(loc*7/8):min(length(dist),floor(loc*9/8))),ceil(loc*7/8));
num8=round(2*durs/quarter);
lengths = cell(length(durs),1);
for n = 1:length(durs)
    switch(num8(n))
        case 1
            lengths(n) = cellstr('eighth note');
        case 2
            lengths(n) = cellstr('quarter note');
        case 3
            lengths(n) = cellstr('dotted quarter note');
        case 4
            lengths(n) = cellstr('half note');
        case 6
            lengths(n) = cellstr('dotted half note');
        case 8
            lengths(n) = cellstr('whole note');
        otherwise
            lengths(n) = cellstr('unknown length');
    end
end
sz=size(notes);
music = cell(sz(1),sz(2)+1);
for r = 1:sz(1)
    for c = 1:sz(2)
        music{r,c}=notes{r,c};
    end
    music{r,end}=lengths{r};
end
tempo = 60*fs/quarter;
end

% helper function to find the average x value in a distribution f.
function x = avex2(f,start)
if (nargin < 2)
    start = 1;
end
x = sum(([1:length(f)]+start-1).*f)/sum(f);
end