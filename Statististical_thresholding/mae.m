function [ e ] = mae( f,g )

e = mean(abs(f - g));

end
