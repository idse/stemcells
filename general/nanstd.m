function y = nanstd(x)

    y = std(x(~isnan(x)));
end