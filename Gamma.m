function output = Gamma(m)
    output = normrnd(0,50) - (m*9.8); %normrand(0,50) gives the zero-mean gaussian random input with rms value of 50
end