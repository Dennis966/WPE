function x_bar = x_bar_gen(half_stfts, len, Lw, f_bin, frame, b, n_mics)
                 
x_bar = zeros(len, 1);
for m = 1:n_mics
    for tau = b:(Lw+b-1)
        if (frame-tau) < 1
            x_bar(Lw*(m-1)+(tau-b+1), 1) = 0;
        else
        x_bar(Lw*(m-1)+(tau-b+1), 1) = half_stfts(m, f_bin, frame-tau);
        end
    end
end

end