%%  %%%%%%%%%%%%%%%%%%%%%%%%% postprocessing of algorithm results %%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%%  %%%%%%%%%%%%%%%%%%%%%%%%% display the ECDF plots of the current run  %%%%%%%%%%%%%%%%%%%%%%%%%
%

foldername = [date '_rotKM_RS_' input.strategy];
    efn = exist(foldername);
    if efn ~= 7
        mkdir(foldername);
    end

for j=1:length(DIMENSION)
dim  = DIMENSION(j);
Eint1(:,j)  = ecdf_data{j}.IntervalFevals./dim;
Etar1(:,j)  = ecdf_data{j}.TargetsPerFevals./(problem.number_of_runs*problem.TarNum);
eval(['Stats(' num2str(j) ',:) = StatsN' num2str(DIMENSION(j)) ';']);
end

p1=figure;
    hndl1=eval(['p' num2str(1)]);
    axes('Parent',hndl1,'Position',[0.15 0.13 0.75 0.70],'LineWidth',1,'FontSize',14);
    set(gcf,'color','w');
    set(gcf,'position',[1 600 800 600]);
    for j = 1:length(DIMENSION)
        dim = DIMENSION(j);
        dis  = 'DisplayName';
        name = ['N' num2str(dim)];
        % for each dimension plot the ratio of reached targets against the
        % needed function evaluations
        eval(['h' num2str(j) '=plot(log(Eint1(:,j))/log(10),Etar1(:,j),dis,name);'])
        hold on
    end
    hleg = legend('show');
    legend boxoff
    title(input.strategy)
    plot([log(problem.budget_multiplier)/log(10) log(problem.budget_multiplier)/log(10)],[0 1],'k','HandleVisibility','off')
    plot([0 ceil(log(problem.budget_multiplier)/log(10))+0.5],...
        [(problem.ConTarNum+1)/problem.TarNum...
        (problem.ConTarNum+1)/problem.TarNum],'--k','HandleVisibility','off')
    xlim([0 ceil(log(problem.budget_multiplier)/log(10))+0.5])
    ylim([0 1])
    xlabel('log(#fevals/dimension)')
    ylabel('target ratio')
    set(hleg,'FontSize',14,'Location','bestoutside');
    grid on

     print(gcf,[foldername '/ecdf_' input.strategy '.pdf'],'-dpdf','-r600','-bestfit')

%%  %%%%%%%%%%%%%%%%%%%%%%%%% median error in the objective space  %%%%%%%%%%%%%%%%%%%%%%%%%
%
%
p2=figure;
    hndl2=eval(['p' num2str(2)]);
    ax2=axes('Parent',hndl2,'Position',[0.15 0.13 0.75 0.70],'LineWidth',1,'FontSize',16);
    set(gcf,'color','w');
    set(gcf,'position',[1 600 800 600]);
    loglog(DIMENSION,abs(Stats(:,6)),'-sb','DisplayName','ES','LineWidth',1,'MarkerSize',8,'LineWidth',1.5)
    hold on
    plot([1 50],[10^-8 10^-8],'--k','DisplayName','final target')
    xlim([1.5 max(DIMENSION)+10])
%     ylim([10^-9 10^-0])
    xlabel('dimension')
    ylabel('median obj. function error')
    grid on
    hleg = legend('show');
    legend boxoff
    set(hleg,'FontSize',16,'Location','NorthWest');
        set(ax2,'FontSize',20,'MinorGridlines','none')
        xticks(DIMENSION)

     print(gcf,[foldername '/medianFitError_' input.strategy '.pdf'],'-r600','-dpdf','-bestfit')

%%  %%%%%%%%%%%%%%%%%%%%%%%%% mean error in the parameter space  %%%%%%%%%%%%%%%%%%%%%%%%%
%
%
p3=figure;
    hndl3=eval(['p' num2str(3)]);
    ax3=axes('Parent',hndl3,'Position',[0.15 0.13 0.75 0.70],'LineWidth',1,'FontSize',16);
    set(gcf,'color','w');
    set(gcf,'position',[1 600 800 600]);
    errorbar(DIMENSION,Stats(1:length(DIMENSION),8),Stats(1:length(DIMENSION),10),':pg','DisplayName', input.strategy,'MarkerSize',8,'LineWidth',1.5)
    hold on
    xlim([1.5 max(DIMENSION)+10])
%     ylim([10^-9 10^-3])
    xlabel('dimension')
    ylabel('mean norm error')
    grid on
    hleg = legend('show');
    legend boxoff
    set(hleg,'FontSize',16,'Location','NorthWest');
    set(ax3,'FontSize',20,'MinorGridlines','none')
    xticks(DIMENSION)
%     yticks([1e-9 1e-8 1e-7 1e-6 1e-5 1e-4 1e-3])
            
     print(gcf,[foldername '/meanNormError_' input.strategy '.pdf'],'-r600','-dpdf','-bestfit')

%%   %%%%%%%%%%%%%%%%%%%%%%%%% print Statistics to TXT-File %%%%%%%%%%%%%%%%%%%%%%%%%
%

    filename = strcat([input.strategy '_RS_on_RotKM_Table.txt']);
    fileID = fopen([foldername '/' filename],'w');        
    fprintf(fileID,'%14s\t %14s\t %14s\t %14s\t %14s\t %14s\t %14s\t %14s\t %14s\n','Dim','Fopt','Fbest','Fmed','\nu','medianErr','FR','meanNormErr','Fevals');
    for j=1:length(DIMENSION)
        fprintf(fileID,'%14.8g\t %14.8g\t %14.8g\t %14.8g\t %14.8g\t %14.8g\t %14.8g\t %14.8g\t %14.8g\n',Stats(j,1:9));
    end
    fclose(fileID);
