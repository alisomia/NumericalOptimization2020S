% set configuration and load the data
clc;
clear;

% set configuration
addpath('./');
addpath(genpath('./database'));
addpath(genpath('./src'));

% load data
load('./database/yale64.mat');
yale64 = fea;
r_yale = 15;
load('./database/orl.mat');
orl = fea;
r_orl = 25;
clear fea gnd
 
