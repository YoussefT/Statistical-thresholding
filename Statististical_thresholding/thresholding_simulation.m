function [mean_error,rate_thresh] = thresholding_simulation(model,A0,thresh_method,alpha,M,j0,Jmax,n,tvec,error_measure,plot_id)
%Perform statistical thresholding on simulated data    
%
%   statistical_thresholding(model,A0,thresh_method,alpha,M,j0,Jmax,n,tvec,error_measure,plot_id)
%
%List of inputs:
%
%   model = 'blocks', 'bumps' or 'trianglesine' intensity model for simulated data.
%   A0 = Value gives a mean intensity over [0,T] equal to 2*A0.
%   thresh_method = 'linear' for no thresholding and projection of the intensity on V_Jmax,
%                   'dm-l' for local hard thresholding (De Miranda and Morettin (2011)) or 
%                   'lrt-l', 'lrt-i' and 'lrt-g' for local, intermediate and global statistical thresholding.
%   alpha = significance level for the statistical thresholding methods, use any value for 'linear' and 'dm-l'.
%   M = number of i.i.d. realizations of the point process per thresholding simulation.
%   j0 = Scale at which the reconstruction of the intensity is preserved (projection on V_j0). 
%        For example j0 = 0 if we want to keep the initial father coefficient only.
%   Jmax = maximum resolution at which the thresholding of the mother coefficients is operated.
%   n = number of simulations to estimate the mean error.
%   tvec = time grid over which the error measure is computed.
%   error = error measure, 'mae' for mean absolute error or 'rmise' for root mean integrated squared error.
%   plot_id = 1 if plot desired, 0 otherwise.


%%%%% trianglesine intensity parameters %%%%%

T = 1; 
nu = 4;
V = 1;
A = 0.2;


%%%%% thresholding parameters %%%%%

hard_coeff = 3;
if strcmp(thresh_method,'linear')
    j_start = Jmax;
else
    j_start = j0;
end
rate_thresh = zeros(1,length(tvec));
error_data = zeros(n,1);


%%%%% thresholding %%%%%

for r=1:n
    ipoints = cell(M,1);
    
    %model selection%
    if strcmp(model,'blocks')
        for m = 1:M
            ipoints{m} = inhomogeneous_poisson1D_blocks(A0);
        end
    elseif strcmp(model,'bumps')
        for m = 1:M
            ipoints{m} = inhomogeneous_poisson1D_bumps(A0);
        end
    elseif strcmp(model,'trianglesine')
        for m = 1:M
            ipoints{m} = inhomogeneous_poisson1D_trianglesine(A0,T,V,nu,A);
        end
    end    
    
    %estimate intensity%
    if strcmp(thresh_method,'linear')
        rate_est = thresholding_linear( ipoints,tvec,T,Jmax );
    elseif strcmp(thresh_method,'lrt-l')
        rate_est = thresholding_fdr( ipoints,tvec,T,j_start,Jmax,alpha );
    elseif strcmp(thresh_method,'lrt-i')
        rate_est = thresholding_rec( ipoints,tvec,T,j_start,Jmax,1-alpha );
    elseif strcmp(thresh_method,'lrt-g')
        rate_est = thresholding_innovation( ipoints,tvec,T,j_start,Jmax,alpha );
    elseif strcmp(thresh_method,'dm-l')
        rate_est = thresholding_hard( ipoints,tvec,T,j_start,Jmax,hard_coeff );
    end 
    rate_thresh = rate_thresh + rate_est;
    
    %compute error measure%
    if strcmp(error_measure,'mae')
        if strcmp(model,'blocks')
            error_data(r) = mae(rate_est,blocks_fun(tvec,A0)');
        elseif strcmp(model,'bumps')
            error_data(r) = mae(rate_est,bumps_fun(tvec,A0)');
        elseif strcmp(model,'trianglesine')
            error_data(r) = mae(rate_est,trianglesine_fun(tvec,A0,T,V,nu,A)');
        end
    elseif strcmp(error_measure,'rmise')
        if strcmp(model,'blocks')
            error_data(r) = rmse(rate_est,blocks_fun(tvec,A0)');
        elseif strcmp(model,'bumps')
            error_data(r) = rmse(rate_est,bumps_fun(tvec,A0)');
        elseif strcmp(model,'trianglesine')
            error_data(r) = rmse(rate_est,trianglesine_fun(tvec,A0,T,V,nu,A)');
        end
    end
    
end 

mean_error = mean(error_data);
rate_thresh = rate_thresh/n ;


%%%%% plots %%%%%

if plot_id
    figure
    plot(tvec,rate_thresh)
    hold on 
    if strcmp(model,'blocks')
        plot(tvec,blocks_fun(tvec,A0))
    elseif strcmp(model,'bumps')
        plot(tvec,bumps_fun(tvec,A0))
    elseif strcmp(model,'trianglesine')
        plot(tvec,trianglesine_fun(tvec,A0,T,V,nu,A))
    end

    if strcmp(thresh_method,'dm-l') || strcmp(thresh_method,'lrt-l') || strcmp(thresh_method,'lrt-i') || strcmp(thresh_method,'lrt-g')
        title(strcat(upper(thresh_method),' thresholding between levels',{' '},num2str(j0),' and',{' '},num2str(Jmax)))
    elseif strcmp(thresh_method,'linear')
        title(strcat('No thresholding and projection on space V_{',num2str(Jmax),'}'))
    end

    legend({'Estimated Intensity','True  Intensity'},'Position',[0.718 0.825 0.1 0.1])
    
    if strcmp(error_measure,'rmise')
        annotation('textbox',[0.3,0,0.4,0.05],...
        'String',strcat('RMISE = ',{' '},num2str(mean_error)))
    elseif strcmp(error_measure,'mae')
        annotation('textbox',[0.3,0,0.4,0.05],...
        'String',strcat('Averaged MAE = ',{' '},num2str(mean_error)))
    end
end

end
