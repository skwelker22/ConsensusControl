%%
%overlapping pairs
%for nnodes = 20
%(2,20)
%(4,18)
%(5,15)
%(8,16)
%for nnodes = 30
%(
%
%
%
%
% n1 = 4; n2 = 14;
ttIx = 1:tt; ttPlot = dT.*(ttIx-1);

nodeLoop=[1:nNodes];
% nodeLoop = [1,2,3,4,5,6,7,8,9,10,15,25];
lenNodeLoop = length(nodeLoop);
uSSS = jet(lenNodeLoop);
hhs = zeros(length(lenNodeLoop),1); nsIx = 1;
legCellTT = cell(length(lenNodeLoop),1);
figure('Name', 'Control Terms');
for ss = nodeLoop
    subplot(231);
    hhs(nsIx) = plot( ttPlot, squeeze(fi_g_plt(1,ss,ttIx)), 'Color', uSSS(nsIx,:));
    hold on;
    subplot(232);
    plot( ttPlot, squeeze(fi_d_plt(1,ss,ttIx)), 'Color', uSSS(nsIx,:)); hold on;
    hold on;
    subplot(233);
    plot( ttPlot, squeeze(fi_gamma_plt(1,ss,ttIx)), 'Color', uSSS(nsIx,:));
    hold on;
    subplot(234);
    plot( ttPlot, squeeze(fi_g_plt(2,ss,ttIx)), 'Color', uSSS(nsIx,:));
    hold on;
    subplot(235);
    plot( ttPlot, squeeze(fi_d_plt(2,ss,ttIx)), 'Color', uSSS(nsIx,:)); hold on;
    hold on;
    subplot(236);
    plot( ttPlot, squeeze(fi_gamma_plt(2,ss,ttIx)), 'Color', uSSS(nsIx,:));
    hold on;
    legCellTT{nsIx} = num2str(ss);
    nsIx = nsIx + 1;
end
subplot(231); lgd = legend(hhs,legCellTT); fontsize(lgd, 9, 'points');
ylabel("f_{g,x}"); 
subplot(232); ylabel("f_{d,x}"); 
subplot(233); ylabel("f_{\gamma,x}");  subplot(234); ylabel("f_{g,y}"); 
subplot(235); ylabel("f_{d,y}");  subplot(236); ylabel("f_{\gamma,y}"); 

%states
figure('Name', 'Position vs. Time');
nsIx = 1;
hhs = zeros(length(lenNodeLoop),1); nsIx = 1;
legCellTT = cell(length(lenNodeLoop),1);
for ii = nodeLoop
    subplot(211);
    hhs(nsIx) = plot( ttPlot, squeeze(xi(1,ii,ttIx)), 'Color', uSSS(nsIx,:)); hold on;
    subplot(212);
    plot( ttPlot, squeeze(xi(2,ii,ttIx)), 'Color', uSSS(nsIx,:)); hold on;
    legCellTT{nsIx} = num2str(ii);
    nsIx = nsIx + 1;
end
%plot navigational point
subplot(211); plot( ttPlot, xd_plot(1,ttIx), 'k' );
subplot(212); plot( ttPlot, xd_plot(2,ttIx), 'k' );
hold off; xlabel('Time [sec]');
subplot(211); ylabel("p_x"); 
lgd = legend(hhs,legCellTT); fontsize(lgd, 9, 'points');
subplot(212); ylabel("p_y");

%states
figure('Name', 'Velocity vs. Time');
nsIx = 1;
hhs = zeros(length(lenNodeLoop),1); nsIx = 1;
legCellTT = cell(length(lenNodeLoop),1);
for ii = nodeLoop
    subplot(211);
    hhs(nsIx) = plot( ttPlot, squeeze(vi(1,ii,ttIx)), 'Color', uSSS(nsIx,:)); hold on;
    subplot(212);
    plot( ttPlot, squeeze(vi(2,ii,ttIx)), 'Color', uSSS(nsIx,:)); hold on;
    legCellTT{nsIx} = num2str(ii);
    nsIx = nsIx + 1;
end
%plot navigational point
subplot(212);
hold off; xlabel('Time [sec]');
subplot(211); ylabel("v_x"); 
lgd = legend(hhs,legCellTT); fontsize(lgd, 9, 'points');
subplot(212); ylabel("v_y");