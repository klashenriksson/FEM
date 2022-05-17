function eoc = EOC(errors, hs)
    eoc = zeros(length(errors)-1, 1);
    for i = 1:length(errors)-1
        eoc(i) = (log(errors(i+1)) - log(errors(i))) / (log(hs(i+1)) - log(hs(i)));
    end
end