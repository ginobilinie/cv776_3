% Author: Johannes L. Schönberger, Jan-Michael Frahm <{jsch, jmf} at cs.unc.edu>

clear;
close all;
clc;

image1 = 'P1180210.mat';
image2 = 'P1180208.mat';

[projMatrix1, projMatrix2, points3D] = reconstructTwoViewModel(image1, image2, true);
