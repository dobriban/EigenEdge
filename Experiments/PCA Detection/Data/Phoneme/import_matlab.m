%% Import data from text file.
% Script for importing data from the following text file:
%   [data path]\phoneme.data
% the data can be downloaded from
% http://statweb.stanford.edu/~tibs/ElemStatLearn/datasets/phoneme.data
% originally it was used in the paper on
%Penalized Discriminant Analysis by Hastie, Buja and Tibshirani. Ann Stat (1995)
%the actual analysis script is in data_analysis.m
%% Initialize variables.
%rename this to the path of the phoneme data
filename = 'C:\Dropbox\Projects\Optimal LSS\Data\Phoneme\phoneme.data';
delimiter = ',';
startRow = 2;

%% Format string:
formatSpec = '%*q%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%s%s%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,'HeaderLines' ,startRow-1, 'ReturnOnError', false);

%% Close the text file.
fclose(fileID);

%% Create output variable
phoneme = [dataArray{1:end-3}];
class = dataArray{end-2};
%% Clear temporary variables
clearvars filename delimiter startRow formatSpec fileID dataArray ans;

%% Save.
save('phoneme.mat','phoneme','class');