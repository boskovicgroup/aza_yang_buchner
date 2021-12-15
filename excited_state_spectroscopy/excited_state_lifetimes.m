%% import data
wavelengths = importdata('spectrum.csv');
delays = importdata('delay.csv');
piperidine = importdata('piperidine.csv');
piptos = importdata('piperidinium_tosylate.csv');
morpholine = importdata('morpholine.csv');
morphtos = importdata('morpholinium_tosylate.csv');
pip_p_OMe_AcPh = importdata('paramethoxy.csv');
pip_p_OMe_AcPh_Tos = importdata('paramethoxy_tosylate.csv');
Cy_AcPh = importdata('acetophenone_cyclohexane.csv');
AcPh = importdata('acetophenone.csv');
%% plot surface plots
surf(wavelengths.data(), delays.data(:,1), piperidine.data(2:end,:));
surf(wavelengths.data(), delays.data(:,2), piptos.data(2:end,:));
surf(wavelengths.data(), delays.data(:,3), morpholine.data(2:end,:));
surf(wavelengths.data(), delays.data(:,4), morphtos.data(2:end,:));
surf(wavelengths.data(), delays.data(:,5), pip_p_OMe_AcPh.data(2:end,:));
surf(wavelengths.data(), delays.data(:,6), pip_p_OMe_AcPh_Tos.data(2:end,:));
surf(wavelengths.data(), delays.data(:,7), AcPh.data(2:end,:));
surf(wavelengths.data(), delays.data(:,8), Cy_AcPh.data(2:end,:));
%% FFT of delayed data
temp1 = [fliplr(piperidine.data(2:end,:)) piperidine.data(2:end,:)].';
piperidine_fft = real(fft(temp1(1:end-1,:)));
start_at_term = 2;
end_at_term = 32;
piperidine_fft = piperidine_fft(start_at_term:end_at_term,:);
[V, E] = eig(piperidine_fft*piperidine_fft.');
[U,~] = eig(piperidine_fft.'*piperidine_fft);
E = fliplr(diag(E).');
g = zeros(size(E));
%compute the indicator function
for i = 1:size(E,2)
    g(i) = E(i)/sum(E(i:end));
end

%% reconstitute
[~,t]=min(g);

reconst = real(ifft(V(:,end-t+2:end)*diag(fliplr(sqrt(E(1:t-1))))*U(:,end-t+2:end).'));
%% Fourier along time
temp2 = [flipud(piperidine.data(2:end,:)); piperidine.data(2:end,:)];
time_fft = real(fft(temp2(1:end-1,:)));
start_at_term = 2;
end_at_term = 32;
time_fft = time_fft(start_at_term:end_at_term,:);
[L, E] = eig(time_fft*time_fft.');
E = fliplr(diag(E).');
g = zeros(size(E));
for i = 1:size(E,2)
    g(i) = E(i)/sum(E(i:end));
end
%% 2D FFT
temp3 = [fliplr([flipud(piperidine.data(2:end,:)); piperidine.data(2:end,:)]) [flipud(piperidine.data(2:end,:)); piperidine.data(2:end,:)]];
twodfft = real(fft2(temp3(1:end-1,1:end-1)));
[L, E] = eig(twodfft*twodfft.');

%% Piperidinyl data at 354.44 nm
subplot(2,2,1)
plot(delays.data(:,2), piptos.data(2:end, find(wavelengths.data() == 354.44)))
hold
plot(delays.data(:,1), piperidine.data(2:end, find(wavelengths.data() == 354.44)))
xlabel('delay [ps]')
ylabel('signal [mOD]')
legend('compound [6H^+] TsO^-','compound 6')
title('354.44 nm slice')
%% morpholine
subplot(2,2,2)
plot(delays.data(:,4), morphtos.data(2:end,find(wavelengths.data() == 348.9500)))
hold
plot(delays.data(:,3), morpholine.data(2:end,find(wavelengths.data() == 348.9500)))
xlabel('delay [ps]')
ylabel('signal [mOD]')
legend('compound 18H^+ TsO^-','compound 18')
title('348.95 nm slice')
%% cyclohexyl and acetophenone
subplot(2,2,3)
plot(delays.data(:,7), AcPh.data(2:end,find(wavelengths.data() == 359.9300)))
hold
plot(delays.data(:,8), Cy_AcPh.data(2:end,find(wavelengths.data() == 359.9300)))
xlabel('delay [ps]')
ylabel('signal [mOD]')
legend('acetophenone','cyclohexylacetophenone, 1')
title('359.93 nm slice')
%% p methoxy
subplot(2,2,4)
plot(delays.data(:,6), pip_p_OMe_AcPh_Tos.data(2:end,find(wavelengths.data() == 355.8100)))
hold
plot(delays.data(:,5), pip_p_OMe_AcPh.data(2:end,find(wavelengths.data() == 355.8100)))
xlabel('delay [ps]')
ylabel('signal [mOD]')
legend('compound [20H^+] OTs^-','compound 20')
title('355.81 nm slice')