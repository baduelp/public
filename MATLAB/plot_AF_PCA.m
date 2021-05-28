function [highscore, score, coeff, explained]=plot_AF_PCA(AFSnorm, popName,popReInd,colorInd, LOI,fig_title, markerInd,markers, samplename, dir, printbool)

% example
% popName={'POP1' 'POP2' 'POP3' 'POP4' 'POP5' 'POP6' 'POP7' 'POP8' 'POP9' 'POP10' 'POP11' 'POP12' 'POP13' 'POP14' 'POP15' 'POP16' 'POP17' 'POP18'};
% plot_AF_PCA(AFSnorm, popName, [1:9 18 10 11 14 15 12 16 17 13],...
%     'name of sample', ...
%     'path\to\Graphs\', 1);

% AFSnorm is a numSITES x numIND matrix of normalized allele frequency per
% site across individuals/populations whose names are given in popName

numSITES = size(AFSnorm,1);
numIND = size(AFSnorm,2);

if isempty(popReInd)
    popReInd=[1:numIND]; % Can be changed to reorder individual colors
end
for i=1:length(popReInd)
popInvInd(i) = find(popReInd==i);
end
popIndNames=popName;%(popReInd);

% PCA
cmap = colormap(parula(max(64, length(popName)*2)));
lcmap=max(64, length(popName)*2);
close(gcf)

if isempty(LOI)
LOI = [1:numIND];% Can be changed to focus on specific individuals
end
% % Center AFS
% AFSmean = mean(AFSnorm(:,LOI), 2);
% AFScentd = AFSnorm(:,LOI)-AFSmean(:)*ones(1,length(LOI));
% no centering of AFS
AFScentd = AFSnorm(:,LOI);
SNPoI=1:numSITES;% Can be changed to focus on subset of sites (for downsampling e.g.)

[coeff,score,latent,tsquared,explained] = pca(AFScentd(SNPoI,:),'NumComponents',2,'Centered',true); %,'algorithm','als'); 

figure1 = figure('Color',[1 1 1]);
set(figure1,'Resize', 'on','PaperOrientation', 'landscape','PaperUnits',...
    'points','PaperSize', [450 300], 'PaperPosition', [0 0 450 300],...
    'InvertHardcopy','off');
hold on
set(gca, 'TickLength',[0 0], 'XColor', 'k', 'YColor', 'k');
axis square
PCind=[1 2];
for pc=1:2
    [~, I]=sortrows(score(:,PCind(pc)));
    highscore(1:20,PCind(pc))=I(end-19:end);
end
h=[];
for i=1:length(LOI)
        h(colorInd(popInvInd(LOI(i))))=plot(coeff(i,PCind(1)), coeff(i,PCind(2)), 'Linestyle', 'none','Marker', markers{markerInd(LOI(i))} ,...
        'MarkerSize', 5, 'MarkerEdgeColor','none', 'MarkerFaceColor',...
        cmap((colorInd(popInvInd(LOI(i)))-1)*floor(lcmap/max(colorInd(popInvInd)))+1,:));
%         text(coeff(i,PCind(1))+2/100*range([min(coeff(:,PCind(1))) max(coeff(:,PCind(1)))]), coeff(i,PCind(2)), popName{colorInd(LOI(i))}, 'FontSize', 8,...
%             'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle')
end
GOI=unique(colorInd(intersect(popReInd,LOI,'stable')),'stable')
hGOI=unique(colorInd(popInvInd(LOI)),'sorted');
legd=legend(h(hGOI),popIndNames(GOI));
set(legd,'EdgeColor',[0 0 0],'Color',[1 1 1],'Location','EastOutside');
box on
xlabel(['PC' num2str(PCind(1)) ' (' num2str(ceil(explained(PCind(1)))) '%)'], 'Fontsize', 10);
ylabel(['PC' num2str(PCind(2)) ' (' num2str(ceil(explained(PCind(2)))) '%)'], 'Fontsize', 10);
if ~isempty(fig_title)
title(fig_title, 'Fontsize', 12, 'FontWeight', 'bold');
end
hold off

if printbool
    currentF = pwd;
    cd(dir)
    print(gcf,  '-painters', '-dpdf','-r600',['PC' num2str(PCind(1)) ' & ' num2str(PCind(2))...
        ' of ' samplename ' centered ' date  '.pdf']);
    cd(currentF)
    close(figure1);
end

