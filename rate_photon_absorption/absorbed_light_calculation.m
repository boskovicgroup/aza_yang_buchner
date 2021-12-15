%% import irradiance data for the lamp; this is given at 0.5 m and 500 W (units mW/m^2*nm)
irradiance = importdata('66142.dat');
% correct to 300 W by multiplying with 3/5
irradiance(:,2) = irradiance(:,2)*3/5;

%% import transmittance data for 3 filters round wavelengths to nearest whole number 
transmittance = importdata('transmittance_filters.csv');
T_data = transmittance.data(); 
T_data(:,2:2:6) = T_data(:,2:2:6)./100;
T_data(:,1:2:5) = round(T_data(:,1:2:5));
T_305_cutoff = flipud(T_data(1:2:end-1,1:2));
T_280_cutoff = flipud(T_data(1:2:end-1,5:6));
T_280_bp = flipud(T_data(1:2:end-1,3:4));

%% truncate irradiance data 200 - 800 nm and correct irradiance for a given filter
filter_low = 200;
filter_high = 800;
irradiance = irradiance(irradiance(:,1)>=filter_low & irradiance(:,1)<=filter_high, :);
corrected_irradiance_280_bp(:,1) = irradiance(:,1);
corrected_irradiance_280_cutoff(:,1) = irradiance(:,1);
corrected_irradiance_305_cutoff(:,1) = irradiance(:,1);
corrected_irradiance_280_bp(:,2) = T_280_bp(:,2).* irradiance(:,2);
corrected_irradiance_280_cutoff(:,2) = T_280_cutoff(:,2).* irradiance(:,2);
corrected_irradiance_305_cutoff(:,2) = T_305_cutoff(:,2).* irradiance(:,2);

%% correct irradiance for distance from the reactor different from 0.5 m for initial data
distance = 0.3;
corrected_irradiance_280_bp(:,2) = corrected_irradiance_280_bp(:,2)*0.5^2/distance^2;
corrected_irradiance_280_cutoff(:,2) = corrected_irradiance_280_cutoff(:,2)*0.5^2/distance^2;
corrected_irradiance_305_cutoff(:,2) = corrected_irradiance_305_cutoff(:,2)*0.5^2/distance^2;

%% input reaction volume in cubic centimeters
reaction_volume = 3; 

% input cross section of the reactor in cm
pathlength = 1;

% calculate the area irradiated by the lamp in m^2
area = reaction_volume/(pathlength*1e4);

% calculate the power per wavelength
power_spectrum_280_bp = corrected_irradiance_280_bp(:,2) * area;
power_spectrum_280_cutoff = corrected_irradiance_280_cutoff(:,2) * area;
power_spectrum_305_cutoff = corrected_irradiance_305_cutoff(:,2) * area;

%% import extinction coefficients for the substrate
% wavelenght in nm versus extinction coefficient 
e = importdata('acetophenone_piperdinium_tosylate_E.xlsx').data();


%% 
% Get wavelegnths and extinctions coeffs
w = e(1:2:end,1);
f = e(1:2:end,2);
% change to W from mW
power_spectrum_280_bp = power_spectrum_280_bp/1000;
power_spectrum_280_cutoff = power_spectrum_280_cutoff/1000;
power_spectrum_305_cutoff = power_spectrum_305_cutoff/1000;

%% enter concentrations of solutions
concs = [2e-4 1.5e-4 1e-4 0.75e-4 0.5e-4];

%% get the integral between power spectrum and extinction coefficients
integ_280_bp = trapz(w*1e-9,power_spectrum_280_bp.*w.*(1-10.^(-f*concs)))';
integ_280_cutoff = trapz(w*1e-9,power_spectrum_280_cutoff.*w.*(1-10.^(-f*concs)))';
integ_305_cutoff = trapz(w*1e-9,power_spectrum_305_cutoff.*w.*(1-10.^(-f*concs)))';

%% Calculate rate of absorbed photons, I_a from three experiments
% V = 3e-6 m^3 volume of the reaction cell
% h = 6.63e-34 Js Planck's constant
% c = 3e8 m/s speed of light
% N_a = 6.022e23 Avogadro's constant

V = 3e-6;
h = 6.63e-34;
c = 3e8;
N_a = 6.022e23;
I_a_280_bp = (1./(1e3*N_a*V*h*c)).*integ_280_bp;
I_a_280_cutoff = (1./(1e3*N_a*V*h*c)).*integ_280_cutoff;
I_a_305_cutoff = (1./(1e3*N_a*V*h*c)).*integ_305_cutoff;


%% Initial rates from three experiments
rates_280_bp = [6.368418244202503e-08 3.6542275202503504e-08 2.9637465265356846e-08 1.1698757124584057e-08 9.004285250063597e-09]'
rates_280_cutoff = [6.19465946065304e-07 4.195134400236391e-07 3.6873950540123195e-07 1.9403939705452139e-07 1.1331072603640417e-07]'
rates_305_cutoff = [1.3907897601118435e-07 9.996745729800316e-08 5.837237780940155e-08 3.703186641711926e-08 2.482543533938912e-08]'


%% Fit straight lines to individual experiments
A1 = [ones(5,1) I_a_280_bp];
params_280_bp = inv(A1'*A1)*A1'*rates_280_bp;
A2 = [ones(5,1) I_a_280_cutoff];
params_280_cutoff = inv(A2'*A2)*A2'*rates_280_cutoff;
A3 = [ones(5,1) I_a_305_cutoff];
params_305_cutoff = inv(A3'*A3)*A3'*rates_305_cutoff;

%% Create figure using function at the end of the script that plots the data and best straight lines
subplot(2,2,1)
plot(irradiance(:,1), log(irradiance(:,2)),'LineWidth',3)
xlabel('wavelength [nm]')
ylabel('ln(irradiance at 0.5 m, 500 W [mW/m^2nm])')

subplot(2,2,2)
plot(T_280_bp(:,1), T_280_bp(:,2),'LineWidth',3)
hold
plot(T_280_cutoff(:,1), T_280_cutoff(:,2),'LineWidth',3)
plot(T_305_cutoff(:,1), T_305_cutoff(:,2),'LineWidth',3)
xlabel('wavelength [nm]')
ylabel('transmittance \times 10^{-2}')
legend('280 \pm 5 bandpass filter', '280-nm high-pass filter', '305-nm high-pass filter')

subplot(2,2,3)
plot(w, f,'LineWidth',3)
xlabel('wavelength')
ylabel('extinction coefficient [M^{-1}cm^{-1}]')

createfigure(I_a_280_bp, rates_280_bp, I_a_280_cutoff, rates_280_cutoff, I_a_305_cutoff, rates_305_cutoff, params_280_bp(2), params_280_bp(1),params_280_cutoff(2), params_280_cutoff(1),params_305_cutoff(2), params_305_cutoff(1))


%% Rates versus powers
rates_power_var = [2.894741396983152e-07 1.907413025019139e-07 1.3907897601118435e-07 1.0728149351079539e-07 4.836971173623953e-08]'
powers = [475 375 300 225 150]'

%%
function createfigure(X1, Y1, X2, Y2, X3, Y3, slope1, int1, slope2, int2, slope3, int3)
%CREATEFIGURE(X1, Y1, X2, Y2, X3, Y3, XData1, YData1, YData2, YData3)
%  X1:  vector of x data
%  Y1:  vector of y data
%  X2:  vector of x data
%  Y2:  vector of y data
%  X3:  vector of x data
%  Y3:  vector of y data
%  XDATA1:  line xdata
%  YDATA1:  line ydata
%  YDATA2:  line ydata
%  YDATA3:  line ydata

%  Auto-generated by MATLAB on 12-Nov-2021 09:54:32


% Create axes

axes1 = subplot(2,2,4);
hold(axes1,'on');

% Create plot

plot(X1,Y1,'DisplayName','280-nm band pass filter',...
    'MarkerFaceColor',[0 0.447058826684952 0.74117648601532],...
    'MarkerSize',10,...
    'Marker','o',...
    'LineStyle','none');

% Create plot
plot(X2,Y2,'DisplayName','280-nm high pass filter',...
    'MarkerFaceColor',[0.850980401039124 0.325490206480026 0.0980392172932625],...
    'MarkerSize',10,...
    'Marker','o',...
    'LineStyle','none');

% Create plot
plot(X3,Y3,'DisplayName','305-nm high pass filter',...
    'MarkerFaceColor',[0.929411768913269 0.694117665290833 0.125490203499794],...
    'MarkerSize',10,...
    'Marker','o',...
    'LineStyle','none');
%range
range = 0:1e-8:7e-7

% Create line
line(range, range*slope1+int1,'Parent',axes1,'MarkerSize',10,'LineWidth',2,...
    'LineStyle','--');

% Create line
line(range, range*slope2+int2,'Parent',axes1,'MarkerSize',10,'LineWidth',2, 'LineStyle','--',...
    'Color',[0.850980401039124 0.325490206480026 0.0980392172932625]);

% Create line
line(range, range*slope3+int3,'Parent',axes1,'MarkerSize',10,'LineWidth',2,'LineStyle','--',...
    'Color',[0.929411768913269 0.694117665290833 0.125490203499794]);

% Create ylabel
ylabel('initial reaction rate [M s^{-1}]');

% Create xlabel
xlabel('rate of photon absorption, I_a [M s^{-1}]');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(axes1,[0 0.7]);
% Uncomment the following line to preserve the Y-limits of the axes
% ylim(axes1,[-6.67716036812981e-08 7.36761765971676e-07]);
% Uncomment the following line to preserve the Z-limits of the axes
% zlim(axes1,[-1 1]);
hold(axes1,'off');
% Set the remaining axes properties
set(axes1,'FontSize',16);
% Create legend
legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.141527138240706 0.733502538071065 0.311501597444089 0.149746192893401],...
    'EdgeColor',[1 1 1]);
end
