% DT2024 - Group 2 - Data arrangement
% Device: HP LP2475w

clear; close all; clc

% Loading required data
load xyz31_1nm              % CIE 1931 2-deg CMFs
load Wavelength_Hamamatsu   % Measured wavelengths

% Selecting relevant CMF data and defining new wavelength range
cmf = xyz31_1nm(21:2:421,2:4); 
lambda_old = wavelength_Hamamatsu;
lambda = 380:2:780;

% Initializing matrix to store spectra data
I = zeros(201, 110);

% First column 
I(:,1) = lambda; 

% Analyzing warm-up data
for i = 1:16
    % Importing warm-up data
    filename = fullfile('data/warm-up/', strcat('ledwarmup', int2str(i), '.txt'));
    data = importdata(filename);
    spectra = interp1(lambda_old, data(:,2), lambda);
    I(:, 1 + i) = spectra;  % Storing spectra data in matrix

    % Calculating XYZ values
    XYZ = calculateXYZ1(spectra, cmf);

    % Storing results
    C{1,1}(i,2) = XYZ(:, 2);    
    C{1,1}(i,1) = (i-1) * 60;  % Storing angle information
end

%% Analyzing uniformity data

for i=1:25
    % Importing uniformity data
    filename = fullfile('data/patches/', strcat('uniformity', int2str(i), '.txt'));
    data = importdata(filename);
    spectra = interp1(lambda_old, data(:,2), lambda);
    I(:, 17 + i) = spectra;  % Storing spectra data in matrix

    % Calculating XYZ and xy values
    XYZ = calculateXYZ1(spectra, cmf);
    xy = XYZtoxy(XYZ);

    % Storing results
    C{1,2}(i,1) = XYZ(:, 2);
    C{1,2}(i,2:3) = xy(:, 1:2);  % Storing xy coordinates
end

%% Analyzing angle data

for i=1:11
    % Importing angle data
    filename = fullfile('data/angle/', strcat('angle', int2str(i), '.txt'));
    data = importdata(filename);
    spectra = interp1(lambda_old, data(:,2), lambda);
    I(:, 42 + i) = spectra;  % Storing spectra data in matrix

    % Calculating XYZ and xy values
    XYZ = calculateXYZ1(spectra, cmf);
    xy = XYZtoxy(XYZ);

    % Storing results
    C{1,3}(i,1) = XYZ(:, 2);
    C{1,3}(i,2:3) = xy(:, 1:2);  % Storing xy coordinates
end

%% Analyzing ramp data

for i=0:17
    % Importing ramp data for red, green, and blue channels
    filename_r = fullfile('data/ramps/', strcat('r', int2str(i), '.txt'));
    filename_g = fullfile('data/ramps/', strcat('g', int2str(i), '.txt'));
    filename_b = fullfile('data/ramps/', strcat('b', int2str(i), '.txt'));

    % Interpolating spectra
    data_r = importdata(filename_r);
    data_g = importdata(filename_g);
    data_b = importdata(filename_b);

    spectra_r = interp1(lambda_old, data_r(:,2), lambda);
    spectra_g = interp1(lambda_old, data_g(:,2), lambda);
    spectra_b = interp1(lambda_old, data_b(:,2), lambda);

    % Storing results
    I(:, 53 + i) = spectra_r;
    I(:, 70 + i) = spectra_g;
    I(:, 87 + i) = spectra_b;

    % Storing spectra data for all channels
    C{1,4}(i+1,:) = spectra_r';
    C{1,4}(i+1+18,:) = spectra_g';
    C{1,4}(i+1+36,:) = spectra_b';
end

%% Analyzing additional ramp data (magenta)

for i=1:6
    % Importing magenta ramp data
    filename_m = fullfile('data/ramps/', strcat('m', int2str(i+2), '.txt'));
    data = importdata(filename_m);
    spectra = interp1(lambda_old, data(:,2), lambda);
    I(:, 104 + i) = spectra;  % Storing spectra data in matrix

    % Storing results
    C{1,5}(i,:) = spectra';  % Storing spectra data for mix colors
end

% Saving data to file
% save("HPLP2475w.mat", "C", "lambda", "cmf", "I")

%% Functions

% Function to calculate Tristimulus values (XYZ)
function [output] = calculateXYZ1(spectra, xyz_curves)
    k = 683;
    output = k * (spectra * xyz_curves);  % [X Y Z]
end

% Function to calculate xy chromaticity coordinates
function [output1] = XYZtoxy(XYZ)
    % Calculate x, y, and z
    x = XYZ(:,1) ./ (XYZ(:,1) + XYZ(:,2) + XYZ(:,3)); 
    y = XYZ(:,2) ./ (XYZ(:,1) + XYZ(:,2) + XYZ(:,3));
    z = XYZ(:,3) ./ (XYZ(:,1) + XYZ(:,2) + XYZ(:,3));

    % Storing results
    output1(:,1) = x;
    output1(:,2) = y;
    output1(:,3) = z;
end
