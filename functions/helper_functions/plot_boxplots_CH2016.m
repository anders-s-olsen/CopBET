function plot_boxplots_CH2016(entropy,tbl,name)
LSDsessions = [1,3];

%%% Figure
figure,hold on,
for LSDses = 1:2
    subplot(1,2,LSDses)
    boxplot(entropy(tbl.session==LSDsessions(LSDses)),tbl.condition(tbl.session==LSDsessions(LSDses))),hold on
    plot([ones(15,1),2*ones(15,1)]',[entropy(strcmp(tbl.condition,'ses-PLCB')&tbl.session==LSDsessions(LSDses)),entropy(strcmp(tbl.condition,'ses-LSD')&tbl.session==LSDsessions(LSDses))]','k-')
    [h,p,ci,stats] = ttest(entropy(strcmp(tbl.condition,'ses-PLCB')&tbl.session==LSDsessions(LSDses)),entropy(strcmp(tbl.condition,'ses-LSD')&tbl.session==LSDsessions(LSDses)));
    title({['LSD dataset (session ',num2str(LSDsessions(LSDses)),'):'],[name,', p=',num2str(p)]})
end