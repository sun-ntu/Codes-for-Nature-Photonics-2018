close all
clear all


%############ enter the file name   #####################
file = ''; % enter the file name in ?

%############ fitting range   ###########################
Tstart = ; % enter the starting and ending time of the calculation (ps)
Tend = ; 

%############ rough echo range    #######################
% enter the time of the each echo
echoStart = ; %ps
echoEnd = ; %ps

echo2Start = ; %ps
echo2End = ; %ps

echo3Start = echo2Start+((echo2Start+echo2End)/2-(echoStart+echoEnd)/2); %ps
echo3End = echo2End+((echo2Start+echo2End)/2-(echoStart+echoEnd)/2); %ps

echo4Start = ; %ps
echo4End = ; %ps
%#########################################################

dat = '.dat';
txt = '.txt';
A = load([file dat]); 
I = A(:,2);
Ioriginal = I;
t = A(:,1);

I = interP(t,I,echoStart,echoEnd);
I = interP(t,I,echo2Start,echo2End);
I = interP(t,I,echo3Start,echo3End);
I = interP(t,I,echo4Start,echo4End);

% ######################### fitting #################################

i = 1;
while t(i) < Tstart
      i = i+1;
end
ts = i;

i = 1;
while t(i) < Tend
      i = i+1;
end
tend = i;


I = I(ts:tend)/(max(I)-min(I)); I1 = I;
t = t(ts:tend);

Ioriginal = Ioriginal(ts:tend)/(max(Ioriginal)-min(Ioriginal));

weighting = ones(numel(I),1);
p = polyfitweighted(t,I,9,weighting);
Iexp = polyval(p,t);


I = I-Iexp;
I2 = Ioriginal-Iexp;

fs = (numel(t)-1)/(max(t)-min(t));
f = fs*(0:(numel(t)-1))/numel(t);
IF = fft(I);

j = numel(f);
while f(j) > 0.6
    IF(j) = 0;
    j = j-1;
end


if mod(numel(f),2)==1
	IF(ceil(numel(f)/2)+1:numel(f)) = rot90(conj(IF(2:ceil(numel(f)/2))),2);
else
	IF(numel(f)/2+2:numel(f)) = rot90(conj(IF(2:numel(f)/2)),2);
end

II = ifft(IF);

Ifinal = I2-II;


com = 'ReBG_';
txt='.txt';

M = zeros(2*numel(t),1);
for i = 1:numel(t)
    M(2*i-1) = t(i);
    M(2*i) = Ifinal(i);
end

fid = fopen([com  file dat],'w');
fprintf(fid,'%7.4f\t%12.8f\r\n',M);
fclose(fid);



figure(1)
plot(t,Ioriginal,t,Iexp,t,I1,'g--')
figure(2)
plot(t,I2,t,II,'g--')
figure(3)
plot(t,Ifinal)