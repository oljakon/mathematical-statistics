function lab1()
    % Выборка объема n из генеральной совокупности Х
    X = [-1.12,-1.06,0.46,-0.39,0.09,-1.44,-1.64,0.86,-0.24,-1.71, ...
        -0.84,-1.19,-0.84,-0.55,-1.11,-1.84,-0.60,-0.92,-0.69,0.23, ...
        0.51,-2.41,-0.53,-1.41,-0.23,-0.89,-0.13,-1.50,0.02,0.27, ...
        -0.75,-0.06,-0.48,0.14,0.20,-2.22,-1.42,-0.54,0.83,-1.77, ...
        -0.10,-0.07,-0.94,-0.13,-1.76,-0.77,-1.26,-0.29,-1.11, ...
        -0.56,1.19,-0.92,-2.02,-1.94,-0.36,-2.09,-2.51,-1.82,0.39, ...
        -2.08,-0.60,-1.38,-1.12,-0.34,0.77,-1.34,0.24,-0.30,-1.67, ...
        -1.50,-0.77,-0.10,-0.39,-0.35,-2.23,-0.84,-0.85,-0.44,-0.20, ...
        -1.76,-0.91,-1.30,-2.03,-2.50,1.08,0.19,0.03,1.17,-0.05, ...
        -2.88,-1.13,-0.05,-1.37,-0.22,0.88,-1.04,-0.52,-1.64,-0.43, ...
        -0.09,-2.44,-0.78,-2.48,-1.16,-0.44,-0.34,-0.60,-0.11,-0.41, ...
        -0.04,-1.09,-1.81,-0.74,-1.07,-1.07,-0.68,-0.36,-0.65,-1.72, ...
        -0.49];
    
    % Минимальное значение
    Mmin = min(X);
    fprintf('Mmin = %f\n', Mmin);
    
    % Максимальное значение
    Mmax = max(X);
    fprintf('Mmax = %f\n', Mmax);
    
    % Размах выборки
    R = Mmax - Mmin;
    fprintf('R = %f\n', R);
    
    % Выборочное среднее
    mu = mean(X);
    fprintf('mu = %f\n', mu);
    
    % Состоятельная оценка дисперсии
    s2 = var(X);
    fprintf('S2 = %f\n', s2);

    % Нахождение количества интервалов
    m = floor(log2(length(X))) + 2;
    
    % Разбиваем выборку на m интервалов от min до max
    [count, edges] = histcounts(X, m, 'BinLimits', [min(X), max(X)]);
    countLen = length(count);
    
    % Интервалы и количество элементов в них
    fprintf('\nИнтервальная группировка значений выборки при m = %d \n', m);
    for i = 1 : (countLen - 1)
        fprintf('[%f : %f) - %d\n', edges(i), edges(i + 1), count(i));
    end
    fprintf('[%f : %f] - %d\n', edges(countLen), edges(countLen + 1), count(countLen));
    
    % Гистограмма
    plotHistogram(X, count, edges, m);
    hold on; 
    % График функции плотности распределения вероятностей нормальной случайной величины
    f(X, mu, s2, m, R);
    figure;
    % График эмпирической функции распределения
    plotEmpiricalF(X);
    hold on;
    % График функции распределения нормальной случайной величины
    F(sort(X), mu, s2, m, R);
    
function plotHistogram(X, count, edges, m)
    h = histogram();
    h.BinEdges = edges;
    h.BinCounts = count / length(X) / ((max(X) - min(X)) / m);
end

function f(X, MX, DX, m, R)
        delta = R/m;
        sigma = sqrt(DX);
        Xn = min(X):delta/20:max(X);
        Y = normpdf(Xn, MX, sigma);
        plot(Xn, Y, 'red');
end

function F(X, MX, DX, m, R)
        delta = R/m;
        Xn = min(X):delta/20:max(X);
        Y = 1/2 * (1 + erf((Xn - MX) / sqrt(2*DX))); 
        plot(Xn, Y, 'red');
end

function plotEmpiricalF(X)  
        [yy, xx] = ecdf(X);
        stairs(xx, yy);
end

end
