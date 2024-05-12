% DT2024 - Group 2 - Data arrangement
% Device: HP LP2475w

clear; close all; clc

load xyz31_1nm
load Wavelength_Hamamatsu

cmf = xyz31_1nm(21:2:421,2:4);              % 201x3 xBar, yBar, zBar
lambda_old = wavelength_Hamamatsu;          % measured wavelengths
lambda = 380:2:780;                         % wavelength 380-780 nm with 2 nm interval
I(:,1) = lambda; 

for i = 1:16
    filename = fullfile('data/warm-up/', strcat('ledwarmup', int2str(i), '.txt'));
    data = importdata(filename);
    spectra = interp1(lambda_old, data(:,2),lambda);
    I(:, 1 + i) = spectra; 

    XYZ = calculateXYZ1(spectra, cmf);

    C{1,1} (i,2) = XYZ(:, 2);    
    C{1,1} (i,1) = (i-1)*60;
end

%%

for i=1:25
    filename = fullfile('data/patches/', strcat('uniformity', int2str(i), '.txt'));
    data = importdata(filename);
    spectra = interp1(lambda_old, data(:,2),lambda);
    I(:, 17 + i) = spectra; 

    XYZ = calculateXYZ1(spectra, cmf);

    xy = XYZtoxy(XYZ);

    C{1,2} (i,1) = XYZ(:, 2);
    C{1,2} (i,2:3) = xy(:, 1:2);
end

%%

for i=1:11
    filename = fullfile('data/angle/', strcat('angle', int2str(i), '.txt'));
    data = importdata(filename);
    spectra = interp1(lambda_old, data(:,2),lambda);
    I(:, 42 + i) = spectra; 

    XYZ = calculateXYZ1(spectra, cmf);

    xy = XYZtoxy(XYZ);

    C{1,3} (i,1) = XYZ(:, 2);
    C{1,3} (i,2:3) = xy(:, 1:2);
end

%%

for i=0:17
    filename = fullfile('data/ramps/', strcat('r', int2str(i), '.txt'));
    data = importdata(filename);
    spectra = interp1(lambda_old, data(:,2),lambda);
    I(:, 53 + i) = spectra;

    C{1,4} (i+1,:) = spectra';
end

for i=0:17
    filename = fullfile('data/ramps/', strcat('g', int2str(i), '.txt'));
    data = importdata(filename);
    spectra = interp1(lambda_old, data(:,2),lambda);
    I(:, 70 + i) = spectra;

    C{1,4} (i+1+18,:) = spectra';
end

for i=0:17
    filename = fullfile('data/ramps/', strcat('b', int2str(i), '.txt'));
    data = importdata(filename);
    spectra = interp1(lambda_old, data(:,2),lambda);
    I(:, 87 + i) = spectra;

    C{1,4} (i+1+36,:) = spectra';
end

%%

for i=1:6
    filename = fullfile('data/ramps/', strcat('m', int2str(i+2), '.txt'));
    data = importdata(filename);
    spectra = interp1(lambda_old, data(:,2),lambda);
    I(:, 104 + i) = spectra;

    C{1,5} (i,:) = spectra';
end

save("HPLP2475w.mat", "C", "lambda", "cmf", "I")

%% Functions

%  Calculate Tristimulus value
function [output] = calculateXYZ1(spectra,xyz_curves)
k = 683;
output = k*(spectra*xyz_curves);  % [X Y Z] as output
end

% %  function to calculate x y z chromaticity coordinates 


function [output1]=XYZtoxy(XYZ)

% calculate x and y and z
x=XYZ(:,1)./(XYZ(:,1)+XYZ(:,2)+XYZ(:,3)); 
y=XYZ(:,2)./(XYZ(:,1)+XYZ(:,2)+XYZ(:,3));
z=XYZ(:,3)./(XYZ(:,1)+XYZ(:,2)+XYZ(:,3));

output1(:,1)=x;
output1(:,2)=y;
output1(:,3)=z;
end