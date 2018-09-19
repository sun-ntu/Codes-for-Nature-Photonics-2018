clear all
close all

% enter the calculated frequency (Fabry-Perot characteristic frequency)
fcal = ; % calculated frequency (THz)

% enter the thickness of the samples
Lthin = ;
Lthick = ;

% enter the file name of the first and without-first echo samples 
fileThin = load('');
fileThick = load('');

file1stThin = load('');
file1stThick = load('');


% ############# FFT ##############
zzz = zeros(5*length(fileThin),1);

I1 = fileThin(:,2); I1 = [I1 ; zzz];
I2 = fileThick(:,2); I2 = [I2 ; zzz];
IP1 = file1stThin(:,2); IP1 = [IP1 ; zzz];
IP2 = file1stThick(:,2); IP2 = [IP2 ; zzz];
t = (0:numel(I1)-1)*(0.1/6);

fs = (numel(t)-1)/(max(t)-min(t));
f = fs*(0:(numel(t)-1))/numel(t);

i=1;
delf = (f(2)-f(1))/2;
while abs(f(i)-fcal) > delf
    i = i+1;
end
fcalInd = i;


IF1 = abs(fft(I1));
IF2 = abs(fft(I2));
IPF1 = abs(fft(IP1));
IPF2 = abs(fft(IP2));


Ithin = IF1;
Ithick = IF2;

Ithin = IF1./IPF1;
Ithick = IF2./IPF2;
% Ithin = Ithin*Ithick(1)/Ithin(1);
% ############ average ###############
flow = fcal - 0.05;
fhigh = fcal + 0.05;

i=1;
delf = (f(2)-f(1))/2;
while abs(f(i)-flow) > delf
    i = i+1;
end
flowInd = i;

i=1;
delf = (f(2)-f(1))/2;
while abs(f(i)-fhigh) > delf
    i = i+1;
end
fhighInd = i;

ratioExp = (sum(Ithin(flowInd:fhighInd)))/(sum(Ithick(flowInd:fhighInd)));

% ############ plots ###############

figure(1)
semilogy(f,IF1,'ro',f,IPF1,'go',f,IF2,'bs',f,IPF2,'ks')
axis([0 1.5 0.001 0.5])
figure(2)
semilogy(f,Ithin,'ro',f,Ithick,'ks')
axis([0 1.5 0 2])


%############setting parameters###############
Vgan = 8020; Dgan = 6150;
Zgan = Vgan*Dgan;

Dsio2 = 2200 ;
Vsio2 = ones(numel(f),1)*5600;
Zsio2 = Vsio2*Dsio2;

r1 = (Zgan-Zsio2)./(Zgan+Zsio2);
t12 = 2*Zsio2./(Zgan+Zsio2);
t21 = 2*Zgan./(Zgan+Zsio2);
r2 = 1;
% ############ solve for attenuation ###############
scatter = 0.6537; % roughness 0.3nm
scatter = 0.3;

ratio = scatter*(-r1(fcalInd));
sol = @(att) (exp(-att*Lthin*2)+ratio*exp(-att*2*Lthin*2)+ratio^2*exp(-att*3*Lthin*2)+ratio^3*exp(-att*4*Lthin*2))...
    /(exp(-att*Lthick*2)+ratio*exp(-att*2*Lthick*2)+ratio^2*exp(-att*3*Lthick*2)+ratio^3*exp(-att*4*Lthick*2));
% sol = @(att) exp(-att*Lthin*2)/exp(-att*Lthick*2);


att = 0;
delsave = 10000;
while att < 1
    del = abs(sol(att)-ratioExp);
%     del = abs(sol(att)-Ithin(fcalInd)/Ithick(fcalInd));
    if delsave > del
        delsave = del;
        attsave = att;
    end
    att = att+0.001;
end
attsave*2

% % ##### output file #####
% % M = zeros(2*numel(tf),1);
% % for i = 1:numel(f)
% %     M(2*i-1) = f(i);
% %     M(2*i) = IF2(i);
% % end
% % 
% % fid = fopen('5s_2345_by1stGaussian_F.txt','w');
% % fprintf(fid,'%7.4f\t%12.8f\r\n',M);
% % fclose(fid);
% % 
% % M = zeros(2*numel(tf),1);
% % for i = 1:numel(f)
% %     M(2*i-1) = f(i);
% %     M(2*i) = IPF2(i);
% % end
% % 
% % fid = fopen('5s_1st_gaussian_F.txt','w');
% % fprintf(fid,'%7.4f\t%12.8f\r\n',M);
% % fclose(fid);