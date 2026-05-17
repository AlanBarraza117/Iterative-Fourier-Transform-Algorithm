% ------------------------------------------------------------------------
% Problem #6: Iterative Fourier Transform Algorithm.
%
% The algorithm is composed of the following steps:
%
% 1.- Apply a FFT from the Diffractive Optical Element (DOE) to output
% plane.
% 2.- Substitute the amplitude of |U_out| with the target amplitude,
% preserving the phase.
% 3.- Apply IFFT to output plane to go back to DOE. 
% 4.- Substitute the amplitude from |U_back| with the original amplitude.
% basically preserve the original amplitude but keep changing the phase.
% 5.- Repeat steps 1 to 4 until the values converge.
%
% Efficiency metrics:
% eta - Diffraction Efficiency: Fraction of power reaching signal window
% UE - Uniformity Error: How much the results differ from the target.
%
% Usage: I only tested my code in matlab not sure if works in Octave.
% ------------------------------------------------------------------------

clear; close all; clc; % I like to always start fresh.

% --------------------------- Parameters ---------------------------------
N = 512; %Size of the grid.
w0_frac = 0.12; %Gaussian Beam Size in proportion to the grid.
spot_sep = round(N/9); %Spot separation in the output plane.
n_iter = 200; %Number of iterations 200 works well.
coords = (0:N-1) - N/2; %Obtaining the center of the grid
[X,Y] = meshgrid(coords, coords);

% -------------------------- Gaussian beam -------------------------------
w0 = w0_frac * N;
A_in = exp( -(X.^2 + Y.^2) / w0^2); %This is a gaussian beam!

% -------------------------- Target Amplitude ----------------------------
% I want to create a grid of 9 Gaussian Beams equally distributed.
r_spot = 4.0; % Size of radius from a spot in pixels.
target = zeros(N,N);
for i = -1:1
    for j = -1:1 % Assigning every position to the spots on the grid.
        cx = i*spot_sep;
        cy = j*spot_sep;
        % Target amplitude this is the output we want to create with IFTA
        target = target + exp(-((X-cx).^2 + (Y-cy).^2) / (2*r_spot^2));
    end
end
target = target/sum(target(:)); % Normalize the quantities.
T_amp = sqrt(target); % Amplitude target.

% -------------------------- Helping tools -------------------------------
% This parameter helps locate pixels with significant energy.
signal_window = target > max(target(:))/10000;
% Customized FFT and IFFT to be in 2D and centered
fft2c = @(u) fftshift( fft2( ifftshift(u) ));
ifft2c = @(U) fftshift( ifft2( ifftshift(U) ));

% ------------------------- IFTA Process ---------------------------------
phi_doe = 2*pi * rand(N,N); % Start with a random phase grid

% This 2 lists will help me track the performance of the algorithm.
% I did use them a lot while writing and testing the program so I will 
% keep them as an extra for the problem.
eta_hist = zeros(1, n_iter); 
ue_hist = zeros(1, n_iter);

fprintf('Starting IFTA for (%d iterations)... \n', n_iter);

for iteration = 1:n_iter

    % Step #1
    U_out = fft2c(A_in .* exp(1j * phi_doe));
    I_out = abs(U_out).^2;

    % Update the metrics:
    I_sig      = I_out(signal_window);
    eta        = sum(I_sig) / sum(I_out(:));  % diffraction efficiency
    T_sig      = target(signal_window);
    scale      = sum(I_sig) / sum(T_sig);
    ue         = sum(abs(I_sig - scale*T_sig)) / (2 * sum(I_sig));
    eta_hist(iteration) = eta;
    ue_hist(iteration)  = ue;

    % Step #2 Substitute the new values for the phase space.
    phase_out = angle(U_out); % Phase extraction
    U_out_c = T_amp .* exp(1j * phase_out); % Keep the original amplitude

    % Step #3 Back-Propagation
    U_back = ifft2c(U_out_c);

    % Step #4 Extract the new phase for the DOE
    phi_doe = angle(U_back);
end


% My performance metrics.
fprintf('\nFinal:  eta = %.1f%%,  UE = %.2f%%\n', ...
        eta_hist(end)*100, ue_hist(end)*100); % <- Transform to percentage

% ------------------ Intensity and DOE Plots -----------------------------
U_final = fft2c( A_in .* exp(1j * phi_doe) );
I_final = abs(U_final).^2;
I_final = I_final/max(I_final(:)); % Normalize so the plot looks good.

I_in_plot    = A_in.^2;
I_in_plot    = I_in_plot / max(I_in_plot(:));
 
target_plot  = target / max(target(:));
 
% Crop central region for output / target display
crop = round(N/4);
sl   = (N/2 - crop + 1):(N/2 + crop);   % MATLAB 1-based slice

% Figure 1:  Input intensity
figure(1);
imagesc(I_in_plot);  colormap(gca,'hot');  axis image;  colorbar;
title('Input intensity |A_{in}|^2  (Gaussian beam)', 'FontSize', 12);
xlabel('x (pixels)');  ylabel('y (pixels)');

% Figure 2:  Target intensity
figure(2);
imagesc(target_plot(sl,sl));  colormap(gca,'hot');  axis image;  colorbar;
title('Target intensity  (3×3 spot array)', 'FontSize', 12);
xlabel('x (pixels)');  ylabel('y (pixels)');

% Figure 3:  DOE phase
figure(3);
imagesc(phi_doe, [-pi, pi]);  colormap(gca,'hsv');  axis image;  colorbar;
title('DOE phase  \phi_{DOE}  (rad)', 'FontSize', 12);
xlabel('x (pixels)');  ylabel('y (pixels)');

% Figure 4:  Reconstructed output intensity
figure(4);
imagesc(I_final(sl,sl));  colormap(gca,'hot');  axis image;  colorbar;
title(sprintf('Reconstructed output intensity  (\\eta=%.1f%%, UE=%.2f%%)', ...
    eta_hist(end)*100, ue_hist(end)*100), 'FontSize', 12);
xlabel('x (pixels)');  ylabel('y (pixels)');

% Figure 5:  Convergence
figure(5);
iters = 1:n_iter;
plot(iters, eta_hist*100, 'b-', 'LineWidth', 2);  hold on;
plot(iters, ue_hist*100,  'r-', 'LineWidth', 2);
xlabel('Iteration');  ylabel('(%)');
legend('Efficiency \eta (%)', 'Uniformity error UE (%)');
title('IFTA convergence', 'FontSize', 12);
grid on;


%-------------------------------------------------------------------------
%
% Code made by:
%               Jose Alan Barraza Villaverde
%               iPSRS (Machine Learning and Photonics)
%
%-------------------------------------------------------------------------