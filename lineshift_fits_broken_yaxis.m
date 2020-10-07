clearvars; clc;

data_mat = csvread('broken_yaxis_plot.csv');
temperature = data_mat(:,1);
position = data_mat(:,2);

figure
BreakPlot(temperature,position,11325,12280,'Line',2000)
