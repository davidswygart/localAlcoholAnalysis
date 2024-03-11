function plotBpod(bpod)







% 
% x = bpod.microInjectionStart;
% plot([x,x], y,'r')
% text(x,y(1), 'inj')
% x = bpod.sipperStart;
% plot([x,x], y, 'r')
% text(x,y(1), 'sip')
% x = bpod.tailStart;
% plot([x,x], y, 'r')
% text(x,y(1), 'tail')
% x = bpod.postInjectionStart;
% plot([x,x], y, 'r')
% text(x,y(1), 'post-inj')
% 
% x = bpod.dropDeployed;
% y = ones(length(x),1) *y(2); 
% scatter(x,y, 'r.')




    %% get bpod times (s) 
    baseline = bpod.baselineStart;
    inj = bpod.microInjectionStart;
    postInj = bpod.postInjectionStart;
    sip = bpod.sipperStart;
    tail = bpod.tailStart;
    tailEnd = tail+ 10*60;
    drop = bpod.dropDeployed;

    %% plot red vertical lines seperating bpod phases
    y = ylim();
    y = [y(2),y(2)*.95];
    plot([baseline,baseline], y,'r')
    plot([inj,inj], y,'r')
    plot([postInj,postInj], y,'r')
    plot([sip,sip], y,'r')
    plot([tail,tail], y,'r')

    %% plot bpod labels
    y = y(1);
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
end