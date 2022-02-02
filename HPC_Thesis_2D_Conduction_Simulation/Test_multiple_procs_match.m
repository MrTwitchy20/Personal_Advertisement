%% Test to ensure both 1 and 36 proc's
%  obtain identical results

opts = detectImportOptions('1_Strong_Scaling_1000x1000_.001__decimal.txt'); %testfile_binary.txt
A = readmatrix('1_Strong_Scaling_1000x1000_.001__decimal.txt',opts); %< Set Matrix to file values

opts = detectImportOptions('36_Strong_Scaling_1000x1000_.001__decimal.txt'); %testfile_binary.txt
B = readmatrix('36_Strong_Scaling_1000x1000_.001__decimal.txt',opts); %< Set Matrix to file values

%opts = detectImportOptions('Numerical_Solution_test_high_toll_decimal.txt'); %testfile_binary.txt
%M = readmatrix('Numerical_Solution_test_high_toll_decimal.txt',opts); %< Set Matrix to file values

figure()
Error = abs(A-B);
contourf(Error, 200, 'linecolor','none');
colormap(jet(256));
title('2D Numerical Error between 1 & 36 procs');
ylabel('y (nodes)');
xlabel('x (nodes)');
zlabel('Error');
colorbar;

