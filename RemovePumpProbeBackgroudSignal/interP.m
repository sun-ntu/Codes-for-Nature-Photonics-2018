function Iout = interP(t,I,ts,te)

% roughly remove the echo from raw data by defining the fitting weighting
tt = abs(t(1)-t(2))/2;

i = 1;
while ts-t(i) > 2*tt
	i = i+1;
end

if abs(ts-t(i)) < tt
	ts = t(i);
	indst1 = i;
else 
	ts = t(i+1);
	indst1 = i+1;
end

i = 1;
while te-t(i) > 2*tt
	i = i+1;
end
if abs(te-t(i)) < tt
	te = t(i);
	indend1 = i;
else 
	te = t(i+1);
	indend1 = i+1;
end
% weighting = ones(numel(I),1);
% weighting(indst:indend) = echoweight;

for i = 1:(indend1-indst1)
	I(i+indst1) = I(indst1)+i*(I(indend1)-I(indst1))/(indend1-indst1);
end

Iout = I;