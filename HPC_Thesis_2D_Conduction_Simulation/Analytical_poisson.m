% \filename Analytical_poisson.m
% \brief    This Mathlab function is used to generate an Analytical 2D
%           plot of temperature distribution. It can be used to read
%           in from .txt files and plot and compare numerical solutions.
% \author   F.OSuibhne
% \date     08.09.21

%% Analytical Solution X,Y Steady State with Boundary Conditions
x_dim = 1; %< Length of x dimension (m)
y_dim = 1; %< Length of y dimension (m)
num_x = 500; %< no. x grid points
num_y = 500; %< no. y grid points
dx = x_dim/(num_x-1); %< x grid spacing 
dy = y_dim/(num_y-1); %< y grid spacing
x = 0:dx:x_dim; % starting at 0 on x axis to 1 divided by dx gaps
y = 0:dy:y_dim;
Grid = zeros(num_y, num_x); %< Initinize grid of 0's
Tm= 100;    % Boudary Conditions
T1= 0;  % Where T1=T2=T3

for i= 1:num_y
    for j= 1:num_x
        Grid(i,j) = Tm*((sinh((pi*y(i))/x_dim)/sinh(pi*y_dim/x_dim)))*sin((pi*x(j))/x_dim) +T1; %< analytical function
    end
end

%% Surface Plot Analytical
figure()
[Y,X] = meshgrid(0:dy:1,0:dx:1);   % Changes plot from node numbers->
%dimension
h= surf(Y,X,Grid,'FaceAlpha',0.9);
set(h,'LineStyle','none')
colorbar;
title('Analytical 2D Steady State');
ylabel('x (m)');
xlabel('y (m)');
zlabel('Temp (°C)');

    
%% Contour Plot Analytical
figure()
%maxGrid = max(Grid,[],'all');
%Norm_An = Grid/maxGrid;
%contourf(Norm_An ,200, 'linecolor','none');
contourf(Y,X,Grid ,200, 'linecolor','none');
colormap(jet(256));
title('Analytical 2D Steady State Temp Distribution');
ylabel('x (m)');
xlabel('y (m)');
zlabel('Temp (°C)');
colorbar;


%% Numerical Solution X,Y Steady State with Boundary Conditions
% .txt is read in with delimiters "spaces" & as Type Cell array

%% Surface Plot Numerical:
figure()
opts = detectImportOptions('Numerical_Sol_500x500_tol.00001.txt'); %testfile_binary.txt
M = readmatrix('Numerical_Sol_500x500_tol.00001.txt',opts); %< Set Matrix to file values
h = surf(Y,X,M,'FaceAlpha',0.9);
set(h,'LineStyle','none')
colorbar;
title('Numerical 2D Steady State Tol:.00001');
ylabel('y (m)');
xlabel('x (m)');
zlabel('Temp (°C)');
%axis([0 500 0 500]);
%% Contour Plot Numerical
figure()
contourf(Y,X,M ,200, 'linecolor','none');
colormap(jet(256));
title('Numerical 2D Steady State Temp Distribution Tol:.00001');
ylabel('y (m)');
xlabel('x (m)');
zlabel('Temp (°C)');
colorbar;
%axis([0 500 0 500]);
%% Error Contour Plot Norms
figure()
maxGrid = max(M,[],'all');
Norm_Num = M/maxGrid;
maxGrid = max(Grid,[],'all');
Norm_An = Grid/maxGrid;
Error = Norm_An-Norm_Num;
contourf(Y,X,Error,200, 'linecolor','none');
colormap(jet(256));
title('Error Between 2D Numerical & Analytical');
ylabel('y (m)');
xlabel('x (m)');
zlabel('Error');
colorbar;
%% Error Contour Plot abs
figure()
Error = abs(Grid-M);
contourf(Y,X,Error,200, 'linecolor','none');
colormap(jet(256));
title('Error 2D Numerical & Analytical Tol=.00001');
ylabel('y (m)');
xlabel('x (m)');
zlabel('Error');
colorbar;