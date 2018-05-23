function [error_data,rate_thresh] = thresholding_with_data(point_data,T,thresh_method,alpha,j0,Jmax,tvec,plot_id)
%Perform statistical thresholding on user input data    
%
%   thresholding_with_data(point_data,thresh_method,alpha,j0,Jmax,tvec,error_measure,plot_id)
%
%List of inputs:
%
%   point_data = point process data input. Multiple realizations of the same process must be regrouped in a cell.               
%   thresh_method = 'linear' for no thresholding and projection of the intensity on V_Jmax,
%                   'dm-l' for local hard thresholding (De Miranda and Morettin (2011)) or 
%                   'lrt-l', 'lrt-i' and 'lrt-g' for local, intermediate and global statistical thresholding.
%   alpha = significance level for the statistical thresholding methods, use any value for 'linear' and 'dm-l'.
%   j0 = Scale at which the reconstruction of the intensity is preserved (projection on V_j0). 
%        For example j0 = 0 if we want to keep the initial father coefficient only.
%   Jmax = maximum resolution at which the thresholding of the mother coefficients is operated.
%   tvec = time grid over which the error measure is computed.
%   plot_id = 1 if plot desired, 0 otherwise.


%%%%% Data preprocessing %%%%%

if not(iscell(point_data))
    point_cell{1} = point_data;
else
    point_cell = reshape(point_data, max(size(point_data)),1);
end
    

%%%%% thresholding parameters %%%%%

hard_coeff = 3;
if strcmp(thresh_method,'linear')
    j_start = Jmax;
else
    j_start = j0;
end
rate_thresh = zeros(1,length(tvec));
error_data = 0;


%%%%% thresholding %%%%%

if strcmp(thresh_method,'linear')
    rate_thresh = thresholding_linear( point_cell,tvec,T,Jmax );
elseif strcmp(thresh_method,'lrt-l')
    rate_thresh = thresholding_fdr( point_cell,tvec,T,j_start,Jmax,alpha );
elseif strcmp(thresh_method,'lrt-i')
    rate_thresh = thresholding_rec( point_cell,tvec,T,j_start,Jmax,1-alpha );
elseif strcmp(thresh_method,'lrt-g')
    rate_thresh = thresholding_innovation( point_cell,tvec,T,j_start,Jmax,alpha );
elseif strcmp(thresh_method,'dm-l')
    rate_thresh = thresholding_hard( point_cell,tvec,T,j_start,Jmax,hard_coeff );
end 


    

%%%%% plots %%%%%

if plot_id
    figure
    plot(tvec,rate_thresh)
    if strcmp(thresh_method,'dm-l') || strcmp(thresh_method,'lrt-l') || strcmp(thresh_method,'lrt-i') || strcmp(thresh_method,'lrt-g')
        title(strcat(upper(thresh_method),' thresholding between levels',{' '},num2str(j0),' and',{' '},num2str(Jmax)))
    elseif strcmp(thresh_method,'linear')
        title(strcat('No thresholding and projection on space V_{',num2str(Jmax),'}'))
    end

    legend({'Estimated Intensity'},'Position',[0.718 0.825 0.1 0.1])
end

end
