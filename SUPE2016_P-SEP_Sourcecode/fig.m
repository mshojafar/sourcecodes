function fig(packets_TO_BS, r, n, yLabel, Title)
% plot data vs. round

    figure(3);
    subplot(1, 3, n);
    plot(1:r, datapackets_TO_BS);
    xlabel('Round');
    ylabel(yLabel);
    title(Title);
end