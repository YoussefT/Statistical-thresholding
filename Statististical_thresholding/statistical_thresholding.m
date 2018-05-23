function [ mean_error,rate_thresh ] = statistical_thresholding( varargin )

data_mode = varargin {1};

if strcmp(data_mode,'data')
    if nargin ~= 9
        error('Missing input argument, see help of thresholding_with_data function');
    else
        if varargin{6}>varargin{7} || varargin{6}<0
            error('j0 must be non-negative and less than or equal to Jmax');
        else
            if ismember(varargin{4},{'lrt-l','lrt-i','lrt-g'}) && check_min_count(varargin{2},varargin{3},varargin{7}) 
                warning('This choice of Jmax might lead to unreliable likelihood ratio test results, M or minimal sum of event counts should be greater than or equal to 50')
                [mean_error,rate_thresh] = thresholding_with_data(varargin{2:end});
            else
                [mean_error,rate_thresh] = thresholding_with_data(varargin{2:end});
            end
        end
    end
elseif strcmp(data_mode,'simulation')
    if nargin ~= 12
        error('Missing input argument, see help of thresholding_simulation function');
    else
        if varargin{7}>varargin{8} || varargin{7}<0
            error('j0 must be non-negative and less than or equal to Jmax');
        else
            if ismember(varargin{4},{'lrt-l','lrt-i','lrt-g'}) && varargin{6}*1.4*varargin{3}/(2^(varargin{8}+1)) < 50 && varargin{6}< 50
                warning('This choice of parameters might lead to unreliable likelihood ratio test results, M or M*1.4*A0/2^(Jmax+1) should be greater than or equal to 50')
                [mean_error,rate_thresh] = thresholding_simulation(varargin{2:end});
            else
                [mean_error,rate_thresh] = thresholding_simulation(varargin{2:end});
            end
        end
    end
else
    error('Mode should be either ''data'' for user input data or ''simulation'' for simulated data');    
end


end
