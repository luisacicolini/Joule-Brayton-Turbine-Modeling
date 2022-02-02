function cycleplot(tablet,tables,air_comb,fuel_comb,run)   
format long
for k=1:5
    namef=run + string(k);
    figure('Name',namef)
    % 0 -> 2t isothermal process
    x=[tables(k,1),tables(k,2)];
    y=[tablet(k,1),tablet(k,2)];
    plot(x,y,'-o');
    hold on
    % 2t -> 3t isentropic process
    x=[tables(k,2),tables(k,3)];
    y=[tablet(k,2),tablet(k,3)];
    plot(x,y,'-o');  
    hold on
    % 3t -> 4t isobaric process
    Tin=tablet(k,3);
    Tfin=tablet(k,4);
    Sin=tables(k,3);
    Sfin=tables(k,4);
    nodes=100;
    deltat=(Tfin-Tin)/nodes;
    deltas=(Sfin-Sin)/nodes;
    xinterval=[Sin:deltas:Sfin Sfin];
    yinterval=[Tin:deltat:Tfin Tfin];
    for st=2:1:nodes
        [cp,cv]=ccalc(yinterval(st-1)+deltat,'polynomial',air_comb(k),fuel_comb(k));
        yinterval(st)=yinterval(st-1)+deltas*(yinterval(st-1)*0.98)/cp; %kinda efficiency
    end
    plot(xinterval,yinterval,'-');
    hold on
    % 4t -> 9t isentropic process
    x=[tables(k,4),tables(k,5)];
    y=[tablet(k,4),tablet(k,5)];
    plot(x,y,'-o');  
    hold on
    title('Turbojet cycle');
    axis([0 2000 0 1500]);    
        %capture the plot as an image
    frame=getframe(gcf);
    im=frame2im(frame);
    [imind, cm]=rgb2ind(im,8);    
    %writing the gif file
    if k==1
        imwrite(imind,cm,'turbojet.gif','gif','LoopCount',inf);
    else 
        imwrite(imind, cm, 'turbojet.gif', 'gif', 'WriteMode','append', 'DelayTime',0.25);
    end
end