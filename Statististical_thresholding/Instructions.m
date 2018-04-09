%% Statistical_thresholding
% Estimation of the intensity from Poisson process data with statistical thresholding methods or
% linear estimation.
%
%% *Syntax*
% For thresholding on user input data:
%%
%
%  statistical_thresholding( data_mode,point_data,T,thresh_method,alpha,j0,Jmax,tvec,plot_id )
%
%%
% For thresholding on simulated data:
%%
%
%  statistical_thresholding( data_mode,model,A0,thresh_method,alpha,M,j0,Jmax,n,tvec,error_measure,plot_id )
%
%% *List of all possible inputs*
% <html>
% <table border=1>
% <tr><td>data_mode</td><td>'data' for user input data</td></tr>
% <tr><td></td><td> 'simulation' for simulated data.</td></tr>
% <tr><td>point_data</td><td> Point process data input. Multiple realizations of the same process must be regrouped in a cell. </td></tr>
% <tr><td> T</td><td>Length of the process specified by user.</td></tr>
% <tr><td>model</td><td>'blocks', 'bumps' or 'trianglesine' intensity model for simulated data.</td></tr>
% <tr><td>A0</td><td>Value gives a mean intensity over [0,T] equal to 2*A0.</td></tr>
% <tr><td>thresh_method</td><td>'linear' for no thresholding and projection of the intensity on V_Jmax.</td></tr>
% <tr><td></td><td>'dm-l' for local hard thresholding (De Miranda and Morettin (2011)).</td></tr>
% <tr><td></td><td>'lrt-l', 'lrt-i' or 'lrt-g' for local, intermediate or global statistical thresholding.</td></tr>
% <tr><td>alpha</td><td>Significance level for the statistical thresholding methods, use any value for 'linear' and 'dm-l'.</td></tr>
% <tr><td>M</td><td>Number of i.i.d. realizations of the point process per thresholding simulation.</td></tr>
% <tr><td>j0</td><td>Scale at which the reconstruction of the intensity is preserved (projection on V_j0). For example j0 = 0 if we want to keep the initial father coefficient only.</td></tr>
% <tr><td>Jmax</td><td>Maximum resolution at which the thresholding of the mother coefficients is operated.</td></tr>
% <tr><td>n</td><td>Number of simulations to estimate the mean intensity and error.</td></tr>
% <tr><td>tvec</td><td>Time grid over which the intensity is estimated and error measure is computed.</td></tr>
% <tr><td>error</td><td>Error measure, 'mae' for mean absolute error or 'rmise' for root mean integrated squared error.</td></tr>
% <tr><td>plot_id</td><td> 1 if plot desired, 0 otherwise.</td></tr></table>
% </html>
%
%% *Author info*
% Created by Youssef Taleb 
%
% Imperial College London
%
% 2018
%