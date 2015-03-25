% Will McFadden (wmcfadden)
% semiflexible polymer cytokinesis simulation

% parameters
d = 120;
L = 1;
Lf = 1;
nf = Lf/L;
f_xl = 0.1;
n_fil =2;
lp = 10;
c = 0.0;
b = 0;

% calculate total segment length from rest length
L0 = L/sqrt(1-L/lp);

% setup positions of endpoints of all segments in all filaments
state = {};
for n=1:n_fil
    x0 = L*floor(rand*d/L);
    dl = Lf;
    if(rand<0.5)
        dl = -Lf;
    end
    x = sort(linspace(x0,x0+dl,nf+1));
    state = {state{:} x};
end
for n=1:length(state)
    q = state{n};
    ind = find(abs(q)<0.001);
    if(~isempty(ind)&&ind~=1)
        if(ind==length(q))
            state(n)=[];
        else
            state{n} = q(ind:end);
        end
        state = {state{:} q(1:ind)+d};
    end
    ind = find(abs(q-d)<0.001);
    if(~isempty(ind)&&ind~=length(q))
        if(ind==1)
            state(n)=[];
        else
            state{n} = q(1:ind);
        end
        state = {state{:} q(ind:end)-d};
    end
end
pl = {};
for n=1:length(state)
    pl = {pl{:} state{n} n*ones(length(state{n}),1)};
end
x = [state{:}]';


% setup filament connection matrix that defines which pairs of xvalues are
% connected in filaments
ind = 0;
ids = x*0;
F = sparse(0,0);
for n=1:length(state)
    nn = length(state{n});
    ids(ind+1:ind+1+nn)=n;
    fadd = sparse([ zeros(nn-1,1) eye(nn-1)]-[  eye(nn-1) zeros(nn-1,1)]);
    F = [F sparse(zeros(size(F,1),nn)); sparse(zeros(nn-1,size(F,2))) fadd];
end

% crosslink filaments by creating a connecting 'C' matrix.  This matrix is
% scaled by "crosslink persistence length" for energy calculation
bins = 0:0.1:d;
C = sparse(0,0);
Cl = [];
for bi=1:length(bins)-1
    inds = find(x>bins(bi)&x<=bins(bi+1));
    for i = 1:length(inds)
        for j = i+1:length(inds)
            if(ids(inds(i))~=ids(inds(j))&&rand<f_xl)
                C = [C; sparse(zeros(1,length(x)))];
                C(end,inds(i))=1;
                C(end,inds(j))=-1;
                Cl = [Cl; x(inds(i))-x(inds(j))];
            end
        end
    end
end
if(isempty(C))
    C =0; 
    Cl = 0;
end

% connect filaments on the edge to resistive spring loads
B = (double(x==0)-double(x==d));
D = d*(x==d);

% energy constant for undisplaced filaments
K = n_fil*nf*(log(L0-L)+L/(L0-L)-1/L)*L0/lp;

x0 = x;

sig = 0.1*L0;
% sample displacements for free energies
% for targ = -0.99:0.01:0.99
x = x0;
oldP = 1;
% fprintf(' %d \n',targ);
errorer = 0.1;
N = 1/errorer^2;
initsteps = 100;
inits = 1000*initsteps;
simsteps = 1000;
totalsteps = simsteps*inits/initsteps;
l = F*x;
stoP = zeros(totalsteps,1);%exp(-(sum(-log(L0-l)-l/(L0-L))*L0/lp + K + sum((C*x - Cl).^2) + b*sum((B.*x+D).^2)));
stol = zeros(totalsteps,1);%mean(l-L)/(L0-L);
stot = zeros(totalsteps,1);%mean(l-L)/(L0-L);
stob = zeros(totalsteps,1);%-mean(B.*x+D);
stoQ = zeros(2,totalsteps);%-mean(B.*x+D);
stox = zeros(length(x),totalsteps);
stol(1) = mean(l);
stoP(1) = 1;
stob(1) = 0;
stoQ(:,1) = [0;0];
stox(:,1) = x;
ind = 2;
l_all = (0:0.01:L0)';
l_sq = repmat(l_all',length(l),1);
% l_all = (0:0.001:L0)';
% l_P = repmat(l_all',length(l_all),1);
% trapz(l_all,trapz(l_all,exp(-0.25*((l_P-l_all(:,ones(1,size(l_P,2))))/sig).^2)'))
% interp_Q = trapz(l_all,exp(-0.5*((l_sq-l_all(:,ones(1,size(l_sq,2))))/sig).^2)');
tic
x0s = zeros(length(x),inits/initsteps);
initind = 1;
fcount = 0;
while initind < inits+1
   fprintf(' %d %d \n',initind, fcount);  
    oldx = x;
    x = x + normrnd(0,sig,size(x));%rand(size(x));
    l = F*x;
    if(sum(l>=L0)==0&&sum(l<=0)==0)
        G = 0;%-sum(log(L0-l)+l/(L0-L)-1./l)*L0/lp + K + c/2*sum((C*x - Cl*exp(-mod(ind,simsteps)/simsteps*10)).^2) + b/2*sum((B.*x+D).^2);
        oldl = F*oldx;
        Q_old = prod(trapz(l_all,exp(-0.25*((l_sq-oldl(:,ones(1,size(l_sq,2))))/sig).^2)'));
        Q_new = prod(trapz(l_all,exp(-0.25*((l_sq-l(:,ones(1,size(l_sq,2))))/sig).^2)'));
        if(rand<exp(-G)/oldP*Q_old/Q_new)
           
            oldP = exp(-G);
        else
            
            x = oldx;
        end
        if(mod(initind,initsteps)==0)
            fprintf(' %d %d \n',initind-1,fcount);  
            x0s(:,(initind)/initsteps)  = x;
            x = x0;
            l = F*x;
            oldP = 1;
        end
        initind=initind+1;
    else
        x = oldx;
        fcount = fcount+1;
%         fprintf(' darn \n'); 
    end
    
    
end
x = x0s(:,1);
whichinit = 2;
while ind < simsteps*inits/initsteps 
    fprintf(' %d %d \n',ind, fcount);  
    oldx = x;
    x = x + normrnd(0,sig,size(x));%rand(size(x));
    l = F*x;
    if(sum(l>=L0)==0&&sum(l<=0)==0)
        G = 0;%-sum(log(L0-l)+l/(L0-L)-1./l)*L0/lp + K + c/2*sum((C*x - Cl*exp(-mod(ind,simsteps)/simsteps*10)).^2) + b/2*sum((B.*x+D).^2);
        oldl = F*oldx;
        Q_old = prod(trapz(l_all,exp(-0.25*((l_sq-oldl(:,ones(1,size(l_sq,2))))/sig).^2)'));
        Q_new = prod(trapz(l_all,exp(-0.25*((l_sq-l(:,ones(1,size(l_sq,2))))/sig).^2)'));
        if(rand<exp(-G)/oldP*Q_old/Q_new)
            stol(ind) = mean(l);
            stoP(ind) = exp(-G);
            stob(ind) = -mean(B.*x+D);
            stot(ind) = mod(ind,simsteps);
            stoQ(:,ind) = [Q_old; Q_new];
            
            stox(:,ind)=x;
            oldP = exp(-G);
        else
            stol(ind) = stol(ind-1);
            stoP(ind) = stoP(ind-1);
            stob(ind) = stob(ind-1);
            stot(ind) = mod(ind,simsteps);
            stoQ(:,ind) = stoQ(:,ind-1);
        
            x = oldx;
            stox(:,ind)=x;
        end
        ind = ind +1;
        if(mod(ind-1,simsteps)==0)
            fprintf(' %d %d !\n',ind-1,fcount);  
            x = x0s(:,whichinit);
            whichinit = whichinit+1;
            l = F*x;
            oldP = 1;
            stol(ind) = mean(l);
            stoP(ind) = 1;
            stob(ind) = -mean(B.*x+D);
            stot(ind) = mod(ind,simsteps);
            stoQ(:,ind) = [0; 0];
            stox(:,ind) = x;
            ind = ind + 1;
            fcount = 0;
        end
    else
        x = oldx;
        fcount = fcount+1;
%         fprintf(' darn \n'); 
    end
    
    
end
toc
% end

stP = stoP;
stl = stol;
stQ = stoQ;
stb = stob;
stx = stox;
stt = stot;

stoP = stP;
stol = stl;
stoQ = stQ;
stob = stb;
stox = stx;
stot = stt;

% divide by partition function to get probabilities
mint = 10;
stoP(stot<mint)=[];
stol(stot<mint)=[];
stoQ(:,stot<mint)=[];
stob(stot<mint)=[];
stox(:,stot<mint)=[];
stot(stot<mint)=[];

r = sortrows([stol stoP]);
l_f = r(:,1);
G_f = r(:,2);
G_f = G_f/trapz(l_f,G_f);

myl = (F*stox)';

l_ff = linspace(min(l_f),max(l_f),100);

[f, ls] = hist(myl(1:mint:end,:),100);
[n,bin] = histc(l_f,l_ff);
G_fff = sparse(1:length(l_f),bin,G_f);
G_ff = full(sum(G_fff)./sum(G_fff~=0));
l_ff(isnan(G_ff))=[];
G_ff(isnan(G_ff))=[];
G_ff = G_ff/trapz(l_ff,G_ff);
figure;
% plot(l_f,G_f,'.');
hold on
% plot(l_ff,G_ff,'r.');
plot(ls,f/trapz(ls,f),'g.');
% reapply updated segment positions to the system state record
ind = 0;
for n=1:length(state)
    state{n}= x(ind+1:ind+length(state{n}));
    ind = ind+length(state{n});
end

% setup for plotting the filaments
pl = {};
for n=1:length(state)
    pl = {pl{:} state{n} n*ones(length(state{n}),1)};
end