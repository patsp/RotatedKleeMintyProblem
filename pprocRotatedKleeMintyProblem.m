clear all;
clc;

foldername = [date '_rotKM_RS_comparison'];
efn = exist(foldername);
if efn ~= 7
  mkdir(foldername);
end

strategies = cell(0, 1);

strategy = struct();
strategy.name = 'RandS';
strategy.plottitle = 'Random Search';
strategy.path = ...
fullfile('.', ...
         [date '_linspaced_rotKM_RS_', ...
          strategy.name]);
strategies{end + 1, 1} = strategy;

set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');

markers = {'+','o','*','h','x','s','d','^','v','>','<','p','.'};
nMarkers = length(markers);

nStrategies = size(strategies, 1);

dimensions = [2, 3];%, 5, 10, 20, 40];
nDimensions = length(dimensions);
for k = 1:nStrategies
  strategies{k, 1}.Eint = cell(nDimensions, 1);
  strategies{k, 1}.Etar = cell(nDimensions, 1);
  strategies{k, 1}.Stats = zeros(nDimensions, 15);
  for j = 1:nDimensions
    dim = dimensions(j);
    load(fullfile(strategies{k, 1}.path, ...
                  ['rotKM_', strategies{k, 1}.name, ...
                   '_Dim', num2str(dim), '.mat']));
    eval(['stats = StatsN', num2str(dim), ';']);
    strategies{k, 1}.Eint{j, 1} = ecdf_data.IntervalFevals ./ dim;
    strategies{k, 1}.Etar{j, 1} = ecdf_data.TargetsPerFevals ./ ...
                                  (1000 * problem.TarNum);
    strategies{k, 1}.Stats(j, :) = stats;
  end
end

problem.budget_multiplier=problem.budget/problem.dim;

for k = 1:nStrategies
  figure('visible', 'off','position',[0, 0, 600, 600]);
  set(gcf,'color','w');
  for j = 1:nDimensions
    dim = dimensions(j);
    dis  = 'DisplayName';
    name = ['N=' num2str(dim)];
    largestX = 6;%ceil(log(problem.budget/dim)/log(10))+1;
    logfevals = log(strategies{k, 1}.Eint{j, 1})/log(10);
    targetRatios = strategies{k, 1}.Etar{j, 1};
    targetRatios(logfevals > largestX) = [];
    logfevals(logfevals > largestX) = [];
    logfevals = [logfevals largestX-0.01 largestX];
    targetRatios = [targetRatios targetRatios(end) targetRatios(end)];
    plot(logfevals,...
         targetRatios, ...
         strcat('-', markers{mod(j-1, nMarkers)+1}), ...
         dis,name, 'LineWidth', 2, ...
         'MarkerSize', 8, ...
         'MarkerIndices', createMarkerIndices(logfevals, 1));
    hold on;
  end
  set(gca, 'FontSize', 12);
  hleg = legend('show');
  legend boxoff;
  title(strategies{k, 1}.plottitle);
  plot([log(problem.budget/dim)/log(10), ...
        log(problem.budget/dim)/log(10)],[0 1],'k','HandleVisibility','off', ...
       'LineWidth', 1.2);
  plot([0 largestX],...
       [(problem.ConTarNum+1)/(problem.ConTarNum+problem.FunTarNum+1)...
                                (problem.ConTarNum+1)/...
        (problem.ConTarNum+problem.FunTarNum+1)],'--k',...
       'HandleVisibility','off', ...
       'LineWidth', 1.2);
  xlim([0 largestX]);
  ylim([0 1]);
  xlabel('log(\#fevals/dimension)');
  ylabel('target ratio');
  set(hleg,'FontSize',12,'Location','southoutside','NumColumns',3);
  grid on;
  set(gcf, 'PaperPositionMode', 'auto');
  pos = get(gcf, 'PaperPosition');
  set(gcf, 'PaperSize', [pos(3) pos(4)]);
  print(gcf, ...
        fullfile(foldername, ['ecdf_' strategies{k, 1}.name '.pdf']), ...
        '-dpdf','-r600','-bestfit');
end

for j = 1:nDimensions
  figure('visible', 'off','position',[0, 0, 600, 600]);
  set(gcf,'color','w');
  dim = dimensions(j);
  for k = 1:nStrategies
    dis1  = 'DisplayName';
    name1 = strategies{k, 1}.plottitle;
    largestX = 6;%ceil(log(problem.budget_multiplier)/log(10))+1;
    len = length(strategies{k, 1}.Etar{j, 1});
    logfevals = log(strategies{k, 1}.Eint{j, 1})/log(10);
    targetRatios = strategies{k, 1}.Etar{j, 1};
    targetRatios(logfevals > largestX) = [];
    logfevals(logfevals > largestX) = [];
    logfevals = [logfevals largestX-0.01 largestX];
    targetRatios = [targetRatios targetRatios(end) targetRatios(end)];
    plot(logfevals, ...
         targetRatios,...
         strcat('-', markers{mod(k-1, nMarkers)+1}), ...
         dis1,name1, ...
         'LineWidth', 2, ...
         'MarkerSize', 8, ...
         'MarkerIndices', createMarkerIndices(logfevals, 1));
    hold on;
  end
  set(gca, 'FontSize', 12);
  plot([log(problem.budget_multiplier)/log(10), ...
        log(problem.budget_multiplier)/log(10)],[0 1], ...
       'k','HandleVisibility','off', 'LineWidth', 1.2);
  title(['N=' num2str(dim)]);
  plot([0 largestX],...
       [(problem.ConTarNum+1)/(problem.ConTarNum+problem.FunTarNum+1)...
                                (problem.ConTarNum+1)/...
        (problem.ConTarNum+problem.FunTarNum+1)],'--k',...
       'HandleVisibility','off', 'LineWidth', 1.2);
  xlim([0 largestX]);
  ylim([0 1]);
  xlabel('log(\#fevals/dimension)');
  ylabel('target ratio');
  grid on;
  hleg = legend('show');
  legend boxoff;
  set(hleg,'FontSize',12,'Location','southoutside','NumColumns',3);
  set(gcf, 'PaperPositionMode', 'auto');
  pos = get(gcf, 'PaperPosition');
  set(gcf, 'PaperSize', [pos(3) pos(4)]);
  print(gcf, ...
        fullfile(foldername, sprintf('ecdf_comparison_dim%04d.pdf', dim)), ...
        '-dpdf','-r600','-bestfit');
end

figure('visible', 'off','position',[0, 0, 600, 600]);
set(gcf,'color','w');
for k = 1:nStrategies
  loglog(dimensions, ...
         max(1e-30, abs(strategies{k, 1}.Stats(:,6))), ...
         strcat('-', markers{mod(k-1, nMarkers)+1}), ...
         'DisplayName',strategies{k, 1}.plottitle,...
         'LineWidth',1,'MarkerSize',8,'LineWidth',2);
  hold on;
end
plot([1 50],[10^-8 10^-8],'--k','DisplayName','final target', ...
    'LineWidth', 1.2);
xlim([1.5 max(dimensions)+10]);
%ylim([10^-9 10^-3]);
yticklabelsQueried = yticklabels();
yticklabelsQueried{1} = '$0$';
yticklabels(yticklabelsQueried);
xlabel('dimension');
ylabel('median obj. function error');
grid on;
hleg = legend('show');
legend boxoff;
set(hleg,'FontSize',12,'Location','southoutside','NumColumns',3);
set(gca,'FontSize',12,'MinorGridlines','none');
xticks(dimensions);
set(gcf, 'PaperPositionMode', 'auto');
pos = get(gcf, 'PaperPosition');
set(gcf, 'PaperSize', [pos(3) pos(4)]);
print(gcf,fullfile(foldername, 'medianFitError_comparison.pdf'), ...
      '-r600','-dpdf','-bestfit');

figure('visible', 'off','position',[0, 0, 600, 600]);
set(gcf,'color','w');
for k = 1:nStrategies
  loglog(dimensions,...
         strategies{k, 1}.Stats(1:length(dimensions),14),...
         strcat('-', markers{mod(k-1, nMarkers)+1}), ...
         'DisplayName', strategies{k, 1}.plottitle, ...
         'MarkerSize',8,'LineWidth',2);
  %{
  errorbar(dimensions,strategies{k, 1}.Stats(1:length(dimensions),14),...
           strategies{k, 1}.Stats(1:length(dimensions),15),...
           strcat('-', markers{mod(k-1, nMarkers)+1}), ...
           'DisplayName', strategies{k, 1}.plottitle, ...
           'MarkerSize',8,'LineWidth',2);
  %}
  hold on;
end
xlim([1.5 max(dimensions)+10]);
%ylim([10^-9 10^-3]);
xlabel('dimension');
ylabel('mean norm error');
grid on;
hleg = legend('show');
legend boxoff;
set(hleg,'FontSize',12,'Location','southoutside','NumColumns',3);
set(gca,'FontSize',12,'MinorGridlines','none');
xticks(dimensions);
%yticks([1e-9 1e-8 1e-7 1e-6 1e-5 1e-4 1e-3]);
set(gcf, 'PaperPositionMode', 'auto');
pos = get(gcf, 'PaperPosition');
set(gcf, 'PaperSize', [pos(3) pos(4)]);
print(gcf,fullfile(foldername, 'meanNormError_comparison.pdf'), ...
      '-r600','-dpdf','-bestfit');

for k = 1:nStrategies
  filename = strcat(strategies{k, 1}.name, '_RS_on_RotKM_Table.txt');
  fileID = fopen(fullfile(foldername, filename),'w');
  fprintf(fileID, ...
          '%14s\t %14s\t %14s\t %14s\t %14s\t %14s\t %14s\t %14s\t %14s\n', ...
          'Dim','Fopt','Fbest','Fmed','\nu', ...
          'medianErr','FR','meanNormErr','Fevals');
  for j=1:length(dimensions)
    fprintf(fileID,['%14.8g\t %14.8g\t %14.8g\t %14.8g\t ', ...
                    '%14.8g\t %14.8g\t %14.8g\t %14.8g\t %14.8g\n'],...
            strategies{k, 1}.Stats(j,[1:7, 14, 10]));
  end
  fclose(fileID);

  filename = strcat(strategies{k, 1}.name, '_RS_on_RotKM_Table.tex');
  fileID = fopen(fullfile(foldername, filename),'w');
  fprintf(fileID, ...
          '%s & %s & %s & %s & %s & %s & %s & %s & %s\\\\\n', ...
          'N','$f_{\text{opt}}$','$f_{\text{best}}$','$f_{\text{med}}$', ...
          '$\nu_{\text{med}}$', '$|f_{\text{med}} - f_{\text{opt}}|$',...
          '$FR$','$||\mathbf{y} - \mathbf{y}_{\text{opt}}||$','meanFevals');
  for j=1:length(dimensions)
    fprintf(fileID,['%.8g & %.8g & %.8g & %.8g & ', ...
                    '%.8g & %.8g & %.8g & %.8g & %.8g\\\\\n'],...
            [strategies{k, 1}.Stats(j,1:5), ...
             abs(strategies{k, 1}.Stats(j,6)),...
             strategies{k, 1}.Stats(j,7), ...
             abs(strategies{k, 1}.Stats(j,14)), ...
             strategies{k, 1}.Stats(j,10)]);
  end
  fclose(fileID);
end

function markerIndices = createMarkerIndices(xs, diff)
  nXs = length(xs);
  markerIndices = [1];
  for k = 2:nXs
    if xs(k) - xs(markerIndices(end)) > diff
      markerIndices = [markerIndices (k - 1)];
    end
  end
end

