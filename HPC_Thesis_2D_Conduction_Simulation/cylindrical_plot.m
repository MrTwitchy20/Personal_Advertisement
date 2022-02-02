%% Polar Plot
% This code was written in order to plot polar, annular grids
% from an rectangular input matrix representing polar space
% the number of nodes in x,y or Thesta,r are used to plot
% The input matrix
% Author    F.OSuibhne 

Num_Theta = 500;    % Theta or y dimension
Num_r = 500; % radial or x dimension
r_i= .02; %internal radius of annulus
r_o= .08; %outer radius of annulus
delta_r= (r_o-r_i)/(Num_r-1);

theta = linspace(0,2*pi,Num_Theta+1);   %< number of nodes
r = linspace(r_i,r_o,Num_r);    %< r_i, r_o number of nodes
[t,r] = meshgrid(theta,r);

opts = detectImportOptions('Annulus_500x500_.00001_Numerical_2_decimal.txt'); %testfile_binary.txt
M = readmatrix('Annulus_500x500_.00001_Numerical_2_decimal.txt',opts); %< Set Matrix to file values

for j= 1:Num_r
    for i= 1:Num_Theta
        Z(j,i) = M(i,j); % inverts data and reads from Numerical source
    end
end

% Program prints MxN however because the solution
% is periodic the below print requires duplication
% of initail column to end on, this allows graph to
% wrap back to beginning

for j=1:Num_r
   Z(j,Num_Theta+1) = M(1,j);  %//set aditional column
end

figure()
[x,y] = pol2cart(theta,r);
h=surf(x,y,Z,'FaceAlpha',0.8);
set(h,'LineStyle','none');
title('Numerical 2D Annulus SS Temp');
ylabel('r (m)');
xlabel('r (m)');
zlabel('Temp (°C)');
colorbar;

figure()
contourf(x,y,Z ,200, 'linecolor','none');
colormap(jet(256));
title('Numerical 2D Annulus SS Temp');
ylabel('r (m)');
xlabel('r (m)');
zlabel('Temp (°C)');
colorbar;
%% Analytical Solution:
% Radial conduction in a tube
% (T-Ti)/(To-Ti)=(ln(r/ri))/(ln(ro/ri));
% T=Ti+(To-Ti)*(ln(r/ri))/(ln(ro/ri));
% Ti=internal, T0=outer temperatures (provided as a boundary condition)
% r_i,r_o are known and r can be a desired grid node

Ti=1;
To=2;
%{
for i = 1:Num_r
        T_Annulus_Ann(i)= Ti+(To-Ti)*((log10((r_i+(i-1)*delta_r)/r_i))/(log10(r_o/r_i)));
end
%}

%% Analytical solution Matrix for 2D 1D heat transfer radially
% 
for j= 1:Num_r
    for i= 1:Num_Theta % Using just the boundary conditions from Numerical 
        T_Annulus_Ann_Matrix(j,i) = Z(1,i) +(Z(Num_r,i) -Z(1,i))*((log10((r_i+(j-1)*delta_r)/r_i))/(log10(r_o/r_i)));
    end
end

for j=1:Num_r
   T_Annulus_Ann_Matrix(j,Num_Theta+1) = T_Annulus_Ann_Matrix(j,1);  %//set aditional column
end

%% Contour Plot Analytical
figure()
contourf(x,y,T_Annulus_Ann_Matrix ,200, 'linecolor','none');
colormap(jet(256));
title('Analytical 2D Annulus SS Temp');
ylabel('r (m)');
xlabel('r (m)');
zlabel('Temp (°C)');
colorbar;

%% Surface Plot Analytical
figure()
h=surf(x,y,T_Annulus_Ann_Matrix,'FaceAlpha',0.8);
set(h,'LineStyle','none')
title('Analytical 2D Annulus SS Temp');
ylabel('r (m)');
xlabel('r (m)');
zlabel('Temp (°C)');
colorbar;

%% Error Contour Plot abs
figure()
Error = abs(Z-T_Annulus_Ann_Matrix);
contourf(x,y,Error,200, 'linecolor','none');
colormap(jet(256));
title('Error 2D Annulus Numerical & Analytical Tol=.00001');
ylabel('r (m)');
xlabel('r (m)');
zlabel('Error');
colorbar;