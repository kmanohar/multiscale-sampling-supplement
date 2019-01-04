%% Setup script to setup path and download NOAA data
% See examples directory to reproduce results in Manohar, Kaiser, Brunton
% and Kutz, "Optimized Sampling for Multiscale Dynamics", https://arxiv.org/abs/1712.05085
%
% Modified 2018-12-31

disp('Adding external/ and its subfolders to matlab path')
addpath(genpath('external'));

disp('Adding src/ and its subfolders to matlab path')
addpath(genpath('src'));

disp('Creating figures subdirectory')
mkdir('figures')

disp('Downloading NOAA SST data from Dropbox URL...')
url = 'https://www.dropbox.com/s/x4511mfaqaca4f4/data.zip?dl=1';
unzip(url,'examples/NOAAtemperature/');

disp('setup external package spgl1')
spgsetup

