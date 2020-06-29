function lab2()
    X = [-1.12,-1.06,0.46,-0.39,0.09,-1.44,-1.64,0.86,-0.24,-1.71,-0.84,...
        -1.19,-0.84,-0.55,-1.11,-1.84,-0.60,-0.92,-0.69,0.23,0.51,-2.41,...
        -0.53,-1.41,-0.23,-0.89,-0.13,-1.50,0.02,0.27,-0.75,-0.06,-0.48,...
        0.14,0.20,-2.22,-1.42,-0.54,0.83,-1.77,-0.10,-0.07,-0.94,-0.13,...
        -1.76,-0.77,-1.26,-0.29,-1.11,-0.56,1.19,-0.92,-2.02,-1.94,-0.36,...
        -2.09,-2.51,-1.82,0.39,-2.08,-0.60,-1.38,-1.12,-0.34,0.77,-1.34,...
        0.24,-0.30,-1.67,-1.50,-0.77,-0.10,-0.39,-0.35,-2.23,-0.84,-0.85,...
        -0.44,-0.20,-1.76,-0.91,-1.30,-2.03,-2.50,1.08,0.19,0.03,1.17,...
        -0.05,-2.88,-1.13,-0.05,-1.37,-0.22,0.88,-1.04,-0.52,-1.64,-0.43,...
        -0.09,-2.44,-0.78,-2.48,-1.16,-0.44,-0.34,-0.60,-0.11,-0.41,...
        -0.04,-1.09,-1.81,-0.74,-1.07,-1.07,-0.68,-0.36,-0.65,-1.72,-0.49];

    % Уровень доверия
    gamma = 0.9;
    
    % Объем выборки 
    n = length(X);
    
    % Точечная оценка матожидания
    mu = mean(X);
    
    % Точечная оценка дисперсии
    sigma2 = var(X);
    
    % Нижняя граница доверительного интервала для матожидания
    muLow = getMuLow(n, mu, sigma2, gamma);
    
    % Верхняя граница доверительного интервала для матожидания
    muHigh = getMuHigh(n, mu, sigma2, gamma);
    
    % Нижняя граница доверительного интервала для дисперсии
    sigma2Low = getSigma2Low(n, sigma2, gamma);
    
    % Верхняя граница доверительного интервала для дисперсии
    sigma2High = getSigma2High(n, sigma2, gamma);
    
    fprintf('mu = %.3f\n', mu);
    fprintf('S2 = %.3f\n', sigma2);
    fprintf('muLow = %.3f\n', muLow);
    fprintf('muHigh = %.3f\n', muHigh);
    fprintf('sigma2Low = %.3f\n', sigma2Low);
    fprintf('sigma2High = %.3f\n', sigma2High);
    
    muArray = zeros(1, n);
    sigma2Array = zeros(1, n);

    muLowArray = zeros(1, n);
    muHighArray = zeros(1, n);
    sigma2LowArray = zeros(1, n);
    sigma2HighArray = zeros(1, n);
    
    for i = 1 : n
        mu = mean(X(1:i));
        sigma2 = var(X(1:i));
        muArray(i) = mu;
        sigma2Array(i) = sigma2;
        muLowArray(i) = getMuLow(i, mu, sigma2, gamma);
        muHighArray(i) = getMuHigh(i, mu, sigma2, gamma);
        sigma2LowArray(i) = getSigma2Low(i, sigma2, gamma);
        sigma2HighArray(i) = getSigma2High(i, sigma2, gamma);
    end
    
    figure 
    hold on;
    plot([1,n], [mu, mu]);
    plot((1:n), muArray);
    plot((1:n), muLowArray);
    plot((1:n), muHighArray);
    xlabel('n');
    ylabel('y');
    legend('mu(x_N)','mu(x_n)','muLow(x_n)','muHigh(x_n)');
    grid on;
    hold off;
    
    figure 
    hold on;
    plot((1:n), (zeros(1, n) + sigma2));
    plot((1:n), sigma2Array);
    plot((1:n), sigma2LowArray);
    plot((1:n), sigma2HighArray);
    xlabel('n');
    ylabel('z');
    legend('sigma(x_N)','sigma(x_n)','sigmaLow(x_n)','sigmaHigh(x_n)');
    grid on;
    hold off;
    
    function muLow = getMuLow(n, mu, s2, gamma)
        muLow = mu - sqrt(s2) * tinv((1 + gamma) / 2, n - 1) / sqrt(n);
    end

    function muHigh = getMuHigh(n, mu, s2, gamma)
        muHigh = mu + sqrt(s2) * tinv((1 + gamma) / 2, n - 1) / sqrt(n);
    end

    function sigma2Low = getSigma2Low(n, s2, gamma)
        sigma2Low = ((n - 1) * s2) / chi2inv((1 + gamma) / 2, n - 1);
    end

    function sigma2High = getSigma2High(n, s2, gamma)
        sigma2High = ((n - 1) * s2) / chi2inv((1 - gamma) / 2, n - 1);
    end

end