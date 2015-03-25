% Will McFadden (wmcfadden)
% simulating relaxation of a collection of nspr springs turning over 
% with time constant tau 

cc = jet(100);
ct = 0;
for tau = [1 0.1]
    nspr = 1000; kappa = 1; sig = 1; totalT = 100*tau;

    figure('Name',['tau = ' num2str(tau)]);
    

    h1= subplot(2,3,1+ct);

    drag = 0.1;
    [t gam] = myturnover(nspr,kappa,tau,drag,sig,totalT);
    plot(t,gam,'color',cc(70,:),'LineWidth',1.5);
    hold on

    drag = 1;
    [t gam] = myturnover(nspr,kappa,tau,drag,sig,totalT);
    plot(t,gam,'color',cc(30,:),'LineWidth',1.5);


    drag = 10;
    [t gam] = myturnover(nspr,kappa,tau,drag,sig,totalT);
    plot(t,gam,'color',cc(50,:),'LineWidth',1.5);

    ylim([0 3]);
    
    h2 = subplot(2,3,2+ct);
    nspr = 100;

    drag = 0.1;
    [t gam] = myturnover(nspr,kappa,tau,drag,sig,totalT);
    plot(t,gam,'color',cc(70,:),'LineWidth',1.5);
    hold on

    drag = 1;
    [t gam] = myturnover(nspr,kappa,tau,drag,sig,totalT);
    plot(t,gam,'color',cc(30,:),'LineWidth',1.5);

    drag = 10;
    [t gam] = myturnover(nspr,kappa,tau,drag,sig,totalT);
    plot(t,gam,'color',cc(50,:),'LineWidth',1.5);

    ylim([0 3]);
    
    h3 = subplot(2,3,3+ct);
    nspr = 50;

    drag = 0.1;
    [t gam] = myturnover(nspr,kappa,tau,drag,sig,totalT);
    plot(t,gam,'color',cc(70,:),'LineWidth',1.5);
    hold on

    drag = 1;
    [t gam] = myturnover(nspr,kappa,tau,drag,sig,totalT);
    plot(t,gam,'color',cc(30,:),'LineWidth',1.5);


    drag = 10;
    [t gam] = myturnover(nspr,kappa,tau,drag,sig,totalT);
    plot(t,gam,'color',cc(50,:),'LineWidth',1.5);

    ylim([0 3]);
ct = 3;
end


% linkaxes([h3,h2,h1])