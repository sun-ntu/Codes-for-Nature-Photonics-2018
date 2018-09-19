clear all
close all

% enter the file name of the first and second echo of thinner sample
thin2 = load('');  % blue (BB)
thick2 = load('');  % black (XX)

% enter the file name of the first and second echo of thicker sample
thin1 = load(''); % BB1
thick1 = load(''); % XX1
noise = 0;

% thick1 --> 256Hz
to = thick1(:,1);
thick1 = thick1(:,2);
for i = 1:numel(to)
    I(2*i-1) = thick1(i);
    if i>2 && i<numel(to)
        I(2*i)= (thick1(i)+thick1(i+1))/2;
    end
end
t = 12:0.05/3:30;
thick1 = [I'; zeros(numel(t)-numel(I),1)];

% thin1 --> 256Hz
to = thick2(:,1);
thick2 = thick2(:,2);
for i = 1:numel(to)
    I(2*i-1) = thick2(i);
    if i>2 && i<numel(to)
        I(2*i)= (thick2(i)+thick2(i+1))/2;
    end
end
t = 12:0.05/3:30;
thick2 = [I'; zeros(numel(t)-numel(I),1)];


N = length(thin2);

BB = zeros(N,1);
XX = zeros(N,1);

BB1 = zeros(N,1);
XX1 = zeros(N,1);

ZZ = zeros(N,1);
for i = 1:N
    t(i) = (0.1/6)*i;
end

for i = 1:length(thin2)
    BB(i) = thin2(i,2);
end

for i = 1:length(thick2)
    XX(i) = thick2(i);
end

for i = 1:length(thin2)
    BB1(i) = thin1(i,2);
end

for i = 1:length(thick2)
    XX1(i) = thick1(i);
end

% for i = 1:floor(length(noise))
%     ZZ(i) = noise(i,2);
% end

XX  = abs(fft(XX))+noise;
BB  = abs(fft(BB))+noise;
XX1  = abs(fft(XX1));
BB1  = abs(fft(BB1));
XXn = XX./XX1;
BBn = BB./BB1;

XXn = XX;
BBn = BB;
XXn = XXn*BBn(1)/XXn(1);

% XX = XX.*ratio;
% XX1 = XX1.*ratio;

% ZZ = abs(fft(ZZ));

D = (4.7)*2; %nm
att = (-log(XXn./BBn)/D)*2; %nm^-1




fs = (numel(t)-1)/(max(t)-min(t));
f = fs*(0:(numel(t)-1))/numel(t);

% fsnoise = (numel(noise(:,1))-1)/(max(noise(:,1))-min(noise(:,1)));
% fnoise = fs*(0:(numel(noise(:,1))-1))/numel(noise(:,1));

% M = zeros(2*numel(f),1);
% for i = 1:numel(f)
%     M(2*i-1) = f(i);
%     M(2*i) = att(i);
% end
% fid = fopen('20s_2nd_F.txt','w');
% fprintf(fid,'%7.4f\t%12.8f\r\n',M);
% fclose(fid);
% 
% M = zeros(2*numel(f),1);
% for i = 1:numel(f)
%     M(2*i-1) = f(i);
%     M(2*i) = BB(i);
% end
% fid = fopen('14s_2nd_F.txt','w');
% fprintf(fid,'%7.4f\t%12.8f\r\n',M);
% fclose(fid);
% 
% M = zeros(2*numel(f),1);
% for i = 1:numel(f)
%     M(2*i-1) = t(i);
%     M(2*i) = thick2(i)*1.16;
% end
% fid = fopen('2nd_20s_normalized.txt','w');
% fprintf(fid,'%7.4f\t%12.8f\r\n',M);
% fclose(fid);



figure(1)
loglog(f,0.1*f.^2,'b-.',f,0.2*f.^4,'k-.',f,att,'ro')
axis([0.1 1.5 0.001 0.2])

% ZZ(1) = 0.001;
% ZZ(numel(ZZ)) = 0.001;

% noisecolor = [0.4 0.4 0.4];



figure(2)
semilogy(f,BB1,'rs','DisplayName','7.6nm 1st echo')
ylabel('FFT power spectrum')
xlabel('Frequency (THz)')
axis([0 2 0.0001 0.3])
hold on
semilogy(f,BB,'ms','DisplayName','7.6nm 2nd echo')
hold on
semilogy(f,XX1,'bv','MarkerFaceColor','b','DisplayName','13.5nm 1st echo')
hold on
semilogy(f,XX,'kv','MarkerFaceColor','k','DisplayName','13.5nm 2nd echo')
hold on
% semilogy(f,ZZ,'g--','MarkerFaceColor','y','DisplayName','noise')
% stem(f,ZZ,':m','fill')
% hold on
% p = fill(f,ZZ,':','facealpha','0.2','DisplayName','noise')
% set(p,'FaceColor',noisecolor)
legend('show')

figure(3)
plot(t+12,thick1*1.16,'k',t+12,thick2*1.16,'b',thin1(:,1),thin1(:,2),'r',thin1(:,1),thin2(:,2),'m')