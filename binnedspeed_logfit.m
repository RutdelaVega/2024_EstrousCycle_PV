%% %%%%%%%% LOG FIT - EXTRACT PARAMETERS AND SAVE FOR STATS %%%%%%%%%% %%

%% 1. Load complete table

filepath = uigetdir();
filename = strcat(filepath, "\BINSallmeandata.xlsx");
bindata = readtable(filename);

%% 2. For all files x remains the same:

binsize = 12;

x = (1:binsize).'; 

%% 3. Loop through bindata:

% 3.1. Vars for loop:
counterstart = 1:binsize:size(bindata, 1)-binsize+1;
counterend = binsize:binsize:size(bindata, 1);
paramstable = nan(size(cat2, 1), 5);

% 3.2. Define model: 

% Fit the data to a logarithmic function y = a + b*log(x)
% Create a model for the fit
model = @(b, x) b(1) + b(2) * log(x);

% Initial guesses for the parameters [a, b]
initialGuess = [1, 1];

a = nan(60, 1); b = a;
R_squared = a; a_pval = a; b_pval = b;

figure(1)
tiledlayout(6, 10)
for ii = 1:size(bindata, 1)/binsize
    tempdata = bindata(counterstart(ii):counterend(ii), 1:5);
    y = tempdata.DFF;
    [params, R, ~, CovB, ~] = nlinfit(x, y, model, initialGuess);
    a(ii) = params(1);
    b(ii) = params(2);
    % Calculate R-squared
    SST = sum((y - mean(y, 'omitnan')).^2); % total sum of squares
    SSE = sum(R.^2, 'omitnan'); % sum of squared errors
    R_squared (ii) = 1 - SSE/SST;
    % Calculate p-values for the parameters
    n = length(y); % number of data points
    p = length(params); % number of parameters
    dof = n - p; % degrees of freedom
    t_stats = params' ./ sqrt(diag(CovB)); % t-statistics for parameters
    pval = 2 * (1 - tcdf(abs(t_stats), dof)); % p-values for parameters
    a_pval(ii) = pval(1); b_pval(ii) = pval(2);

    % Create and save plots
    nexttile
    % Generate fitted curve data
    x_fit = linspace(min(x), max(x), 100); % create a smooth line for the fit
    y_fit = a(ii) + b(ii) * log(x_fit);
    hold on;
    scatter(x, y, 'ro', 'filled'); % plot the original data points
    plot(x_fit, y_fit, 'b-', 'LineWidth', 2); % plot the fitted curve
    hold off;
    titlename = join(string(table2cell(cat2(ii, :))));
    title(titlename);
%     legend('Data', 'Fitted Curve');
    hold off
    clear tempdata y params R CovB  n p dof t_stats SSE SST x_fit y_fit titlename pval
end


paramstable(:, 1) = a; paramstable(:, 2) = a_pval;
paramstable(:, 3) = b; paramstable(:, 4) = b_pval;
paramstable(:, 5) = R_squared; 
 
paramstable = array2table(paramstable, 'VariableNames', ["A (scale)", "A (pval)", "B (slope)", "B (pval)", "R^2"]);

alldata = [cat2, paramstable]; % Cat2 es mi matriz con los ID, puedes guardar directamente paramstable. 

% Se va a guardar en la misma carpeta donde está el archivo de bines
% (filepath)
dataname = strcat(filepath, "\logfit_6cm12bin.xlsx"); % cambiar nombre (si quieres, por ejemplo, logfitPVcreD3.xlsx
% Sustituir alldata por paramstable y se guarda (sin los ID) 
writetable(alldata, dataname)
