clear;
clf;
addpath('../../../');


% open distance
dist = readmatrix('distance.dat'); % meters
radius = readmatrix('radius.dat'); % meters
time = readmatrix('time.dat'); % seconds

% open electric field
fid = fopen('E000.dat', 'r');
A = fread(fid, 'double');
fclose(fid);
complex_A = complex(ones(round(size(A)/2)), ones(round(size(A)/2)));
% make the complex array
Aodd = A(1:2:end);
Aeven = A(2:2:end);
for i = 1:(length(A)/2)
	complex_A(i) = Aodd(i) + Aeven(i)*1i;
end
complex_A = reshape(complex_A, [length(time), length(radius)]);

% size(complex_A)
Intensity = abs(complex_A).^2;
imagesc(Intensity);


wavelength = 800e-9;

% convert lambda to angular frequency
% c = lambda * f
c = 2.998e8; % speed of light
f0 = c / wavelength;
w0 = 2*pi * f0;


% plotting
imagesc(radius*1e3,time*1e15,Intensity);

radius = reshape(radius, [1, length(radius)]);

xv = cat(2, -1*fliplr(radius(1:end)), 0, radius(1:end));
yv = cat(2, -1*fliplr(radius(1:end)), 0, radius(1:end));

% radius_derivative = radius(2:end) - radius(1:end-1);
% figure(1);
% plot(radius_derivative);
% figure(2);
% plot(radius);
% return;

% time % seconds
t_cmc = w0 * time; % periods
t_cmc = reshape(t_cmc, [1, length(t_cmc)]);
% dist

zv = reshape(dist, [1, length(dist)]);

% % % write variables to axes file
% xv = -0.010:0.0020:0.010; %   20um
% yv = -0.010:0.0020:0.010; % x 20um
% zv = -0.025:0.005:0.025;  % x 50um
% t_cmc = 2*pi*(0 : 0.01 : 0.99); % periods
% t_cmc = 2*pi*(-3.0 : 0.01 : 3.0); % periods

% save the axes
save('driving_field/axes', 'xv', 'yv', 'zv', 't_cmc')

% write the field for each point
for z_i=1:length(zv)

	% open electric field
	filename = strcat('E', strcat(num2str(z_i-1, '%03.f'), '.dat'));
	fid = fopen(filename, 'r');
	A = fread(fid, 'double');
	fclose(fid);

	% make the complex array
	complex_A = complex(ones(round(size(A)/2)), ones(round(size(A)/2)));
	Aodd = A(1:2:end);
	Aeven = A(2:2:end);
	for i = 1:(length(A)/2)
		complex_A(i) = Aodd(i) + Aeven(i)*1i;
	end
	complex_A = reshape(complex_A, [length(time), length(radius)]);

	% convert radial coordinates to x y z coordinates
	driving_field = zeros(length(xv),length(yv),1,length(t_cmc));
	for t_i=1:length(t_cmc)
		for xv_i=1:length(xv)
			for yv_i=1:length(yv)
				this_radius = sqrt(xv(xv_i)^2 + yv(yv_i)^2);
				if this_radius < max(radius)
					[m, I] = min(abs(this_radius - radius));
					driving_field(xv_i, yv_i, 1, t_i) = complex_A(t_i, I);
				end
			end
		end
		"time iteration"
		t_i
	end
	"finished the loops"

	% save the driving filed in x, y ,z coordinates
	% XI: x position
	% YI: y position
	% TI: time
	% driving field(XI, YI, C, TI)
	filename = "data_"+string(z_i);
	save('driving_field/'+filename, 'driving_field', '-v7.3');
end

hhgmax = hhgmax_load();
% initialize config struct
config = struct();

config.symmetry = 'rotational';
config.periodic = 0;
config.precomputed_driving_field = 'driving_field';
% config.wavelength = 1e-3; % mm
config.wavelength = wavelength;
config.tau_interval_length = 1; % driving field periods
config.tau_window_length = 0.5; % driving field periods
config.t_window_length = 5; % driving field periods
config.ionization_potential = 12.13; % eV (for Xe)
 propagation_config = struct();
 propagation_config.pressure = 0.1; % bar

 load('xenon_abs.dat');
 propagation_config.transmission_energy = xenon_abs(:,1);
     % photon energy in eV
 propagation_config.transmission = xenon_abs(:,2);
     % for 30 torr, 1 cm

% call harmonic_propagation module which will use the dipole_response module
% for calculating dipole spectra
[z_max, omega, U] = hhgmax.harmonic_propagation(t_cmc, xv, yv, zv, config, propagation_config);

% select one harmonic frequency for plotting
omega_i = 22; % corresponds to 21st harmonic
omega = omega(omega_i)
U = squeeze(U(:,:,1,omega_i));

% plot electric field intensity of harmonic radiation
intensity = abs(U).^2;

colormap(fliplr(hot));
imagesc(xv,yv, intensity / max(max(intensity)));
title(['electric field intensity at z=' num2str(z_max) 'mm'])
xlabel('x [mm]')
ylabel('y [mm]')
