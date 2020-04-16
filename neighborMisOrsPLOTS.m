minAng = 2;
maxAng = 10;

% Plots
misPhases = fieldnames(misXTAL);

% loop to make plots for each phase separately
for k = 1:length(misPhases)
    
    % specified phase
    p = misPhases{k};
    
    % misorientations in crystal coordinates
    mX = misXTAL.(p);
    
    % misorientation axes in specimen coordinates
    mS = misSPEC.(p);
    
    % plotting condition (which angle ranfe to plot?)
    cond = mX.angle<=maxAng*degree&mX.angle>=minAng*degree;
    % PLOTS
    % crystal coordinates
    figure,
    plotAxisDistribution(mX(cond),'antipodal','lower','smooth','halfwidth',10*degree,'figSize','small')
    cb = mtexColorbar('Title','M.U.D.');
    setColorRange([1 max(cb.Limits)]);
    saveas(gcf,sprintf('%s_misO_axes_XTAL.png',p));
    
    % specimen coordinates
    figure,
    plotAxisDistribution(mS(cond),'antipodal','lower','smooth','halfwidth',10*degree,'figSize','small')
    cb = mtexColorbar('Title','M.U.D.');
    setColorRange([1 max(cb.Limits)]);
    saveas(gcf,sprintf('%s_misO_axes_SPEC.png',p));
    
end