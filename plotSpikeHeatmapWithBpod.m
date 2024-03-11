function plotSpikeHeatmapWithBpod(spk, binEdges, bpod, cLabel, cRange)
    imagesc(spk)
    c = colorbar;
    c.Label.String = cLabel;
%     c.Label.Rotation = -90;
%     oldPosition = c.Label.Position;
%     oldPosition(1) = oldPosition(1)+0.9;
%     c.Label.Position = oldPosition;
    set(gca,'box','off')
    ylabel('cluster #')

    clim(cRange)

    oldYlim = ylim();
    labelY = 0.5 - diff(oldYlim)*.14;

    hold on

    %% get bpod times (s) and convert to units of bin time
    binWidth = binEdges(2) - binEdges(1);
    baseline = (bpod.baselineStart - binEdges(1))/ binWidth;
    inj = (bpod.microInjectionStart - binEdges(1))/ binWidth;
    postInj = (bpod.postInjectionStart - binEdges(1)) / binWidth;
    sip = (bpod.sipperStart - binEdges(1))/ binWidth ;
    tail = (bpod.tailStart - binEdges(1)) / binWidth;
    tailEnd = tail+ (10*60/binWidth);
    drop = (bpod.dropDeployed - binEdges(1)) / binWidth;


    %% plot red vertical lines seperating bpod phases
    y=[labelY,0.5];
    plot([baseline,baseline], y,'r')
    plot([inj,inj], y,'r')
    plot([postInj,postInj], y,'r')
    plot([sip,sip], y,'r')
    plot([tail,tail], y,'r')

    %% plot bpod labels
    y = labelY - (labelY-0.5)/2 ;
    x = (inj - baseline)/2 + baseline;
    text(x,y, 'Baseline',HorizontalAlignment='center')

    x = (postInj - inj)/2 + inj;
    text(x,y, 'Inj.',HorizontalAlignment='center')

    x = (sip - postInj)/2 + postInj;
    text(x,y, 'post-Inj',HorizontalAlignment='center')

    x = (tail - sip)/2 + sip;
    text(x,y, 'Sipper',HorizontalAlignment='center')

    x = (tailEnd-tail)/2 + tail;
    text(x,y, 'post-Sip',HorizontalAlignment='center')

    %% plot valve openings
    y = zeros(length(drop),1) + .4;
    scatter(drop,y, 'r.')

    %%
    hold off

    
    ylim([labelY,oldYlim(2)])

end