k = 10;
lambda = [-2*sqrt(k-1):0.1:2*sqrt(k-1)];
plot(lambda, SpectralDensity(lambda,k), '-r', 'linewidth', 2);
ylabel('p(\lambda)');
xlabel('\lambda');
title(['Spectral density, k=',num2str(k)]);
saveas(gcf,['/Users/yangchenye/Downloads/HW3_P1_e_',num2str(k),'.png'])