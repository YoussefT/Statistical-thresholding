function [ e ] = rmse( f,g )

e = sqrt(mean((f - g).^2));

end

