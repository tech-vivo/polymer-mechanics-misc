% Will McFadden (wmcfadden)
function [t stog] = myturnover(nspr,kappa,tau,drag,sig,totalT)
    G0 = kappa;%*(nspr/100)^(11/5);
    spr = zeros(nspr,1);
    totalsteps = totalT/tau;
    gam = 0;
    stog = [];
    stos = [];
    for i = 1:totalsteps
        gamT = (sig/G0 + mean(spr)- gam)*(1-exp(-G0*tau/drag));
        gam = gam + gamT;
        spr(randi(length(spr))) = gam;
        stog = [stog gam];
    end
    t = tau*(1:totalsteps);

end