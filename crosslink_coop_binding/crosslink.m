% Will McFadden (wmcfadden)
% exploring relationship between binding rates and avidity when multiple
% crosslinks bind at a given filament/bundle pair

kon = 10;           % binding rate of individual crosslinker
koff = 0.1;         % unbinding rate
n = 1;              % number of filaments in a bundle
s = 1/360;          %
d = 1/36;
numtry = 10;
cc = jet(numtry);
ind = 1;
for g0 = logspace(2,10,numtry)
    stoc2 = [];
    stoc0 = [];
    stocx = [];
    stonc = [];
    for c00 = logspace(-10,-2,20);
        [t, c] = ode45(@cross_fun, 0:0.1:100,[0 0],[], n,g0, kon, koff, s,g0, c00);
        c2 = c(end,end);
        cx = c2*d/2/g0*n/((n-1) + d/2/g0*n);
        stoc0 = [stoc0 c00];
        stoc2 = [stoc2 c2];
        stocx = [stocx cx];
        stonc = [stonc cx*g0^2];
    end
    plot(log(stoc0), log(stocx),'color',cc(ind,:));
    hold on
    ind = ind+1;
end

