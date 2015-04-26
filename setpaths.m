function [] = setpaths
% 
% This function is used to add the current directory and all of its 
% subdirectories to MATLAB search path in order to allow main functions to use 
% auxiliary procedures stored there.
% =========================================================================
%               PLEASE, RUN THIS FUNCTION FIRST 
%     EACH TIME YOU START A NEW SESSION WITH THE EigenEdge Package
% =========================================================================

addpath(genpath(['.' filesep]));

