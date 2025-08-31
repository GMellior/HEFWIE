% Refine time grid (optional to get smoother dynamics)
i1           = find(Time>=70, 1,'first');
i2           = find(Time>=200,1,'first');
t_start      = Time(i1);
t_end        = Time(i2);
n_old        = i2-i1+1;
n_new        = 2*n_old;                           % Multiplier
Time_refined = linspace(t_start,t_end,n_new + 1)';
Time_refined = Time_refined(2:end);               % Drop first point to avoid duplicate with Time(i1)
Time_new     = [Time(1:i1); Time_refined; Time(i2+1:end)];
dtsvec       = diff(Time_new);
dtsvec       = [dtsvec(1); dtsvec];
T            = Time_new(end);
perds        = length(Time_new);
Time         = Time_new;
