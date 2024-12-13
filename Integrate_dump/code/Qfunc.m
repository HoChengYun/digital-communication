function Q = Qfunc(x)
    % 定義積分步長和積分區間
    dt = 0.001;
    t = x:dt:10; % 積分從 x 到一個較大的上限 (如10)

    % 被積函數的值
    integrand = exp(-t.^2 / 2) / sqrt(2 * pi);

    % 使用梯形法進行數值積分
    Q = sum(integrand) * dt;
end