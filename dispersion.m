clear all 
close all

% enter the file names of the first and second echoes
I1 = load('');
I2 = load('');

fileName = '14s';

% enter the estimated thickness
L = 9;

% t = I1(:,1);
t = 12:0.05/3:100;
fs = (numel(t)-1)/(max(t)-min(t));
f = fs*(0:(numel(t)-1))/numel(t);

I1 = I1(:,2);
I2 = I2(:,2);
I1 = [I1; zeros((numel(t)-length(I1)),1)];
I2 = [I2; zeros((numel(t)-length(I2)),1)];

F1 = fft(I1);
F2 = fft(I2);

angle1 = angle(F1);
angle1 = unwrap(angle1);
angle2 = angle(F2);
angle2 = unwrap(angle2);

angle1 = angle1-angle1(1);
angle2 = angle2-angle2(1);

if angle1(2)>angle2(2)
    v = 2*pi*f'*2*L./((angle1-angle2))*10^3;
else
    v = 2*pi*f'*2*L./((angle2-angle1))*10^3;
end
v(1) = v(2);

L = 5950/v(1)*L;
v = v*5950/v(1);


plot(f,v)
axis([0 1 5200 6200])


% ##### output file #####
i = 1;
while f(i)<1
    i = i+1;
end

M = zeros(i,1);
for j = 1:i
    M(2*j-1) = f(j);
    M(2*j) = v(j);
end

fid = fopen([fileName '_' num2str(L) 'nm.txt'],'w');
fprintf(fid,'%7.4f\t%12.8f\r\n',M);
fclose(fid);