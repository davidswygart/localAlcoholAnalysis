%% Import data from text file
% Script for importing data from the following text file:
%
%    filename: C:\Users\david\Downloads\cluster_info.tsv
%
% Auto-generated by MATLAB on 09-Jan-2024 11:00:01

%% Setup the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 11);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = "\t";

% Specify column names and types
opts.VariableNames = ["cluster_id", "Amplitude", "ContamPct", "KSLabel", "amp", "ch", "depth", "fr", "group", "n_spikes", "sh"];
opts.VariableTypes = ["double", "double", "double", "categorical", "double", "double", "double", "double", "categorical", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, ["KSLabel", "group"], "EmptyFieldRule", "auto");

% Import the data
clusterinfo = readtable("C:\Users\dis006\Downloads\cluster_info.tsv", opts);


%% Clear temporary variables
clear opts



%% count 
summary(clusterinfo.group)
