%Amirhosein_Zaboli_9926843
%Soroush_Rezaei_9926563
%Amir_Jahangard_9925293

%Project1

clc;
clear;
close all; 

branch_and_bus_data
Ybus

method = input('Enter 1 for Gauss-Seidel, 2 for Newton-Raphson, 3 for Fast Decoupled: ');
while method ~= 1 && method ~= 2 && method ~= 3
    fprintf('Invalid Input\n');
    method = input('Enter 1 for Gauss-Seidel, 2 for Newton-Raphson, 3 for Fast Decoupled: ');
end

if method == 1
    maingauss
elseif method == 2
    mainnewton
else 
    maindecouple
end