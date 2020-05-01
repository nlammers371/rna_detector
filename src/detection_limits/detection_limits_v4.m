clear
close all

FigPath = '../../fig/detection_limits_v4/';
mkdir(FigPath);

mol = 6.022e23;
nanomol = 1e-9*6.022e23;
reaction_vol = 2e-5;

wc = 1e-10;
k = 1 / nanomol / reaction_vol; % on rate (roughly diffusion limited)
rcg = 1; % off rate for specific interactions
kc = 200; % cleavage rate

RNP1 = 20*nanomol*reaction_vol;
RNP2 = 0.2*nanomol*reaction_vol;
S0 = 200*nanomol*reaction_vol;
B0 = 2000*nanomol*reaction_vol;


% make figure of "vanilla" NCR 
bc_vec = logspace(-10,0);
a_lim_vec_vanilla = NaN(1,numel(bc_vec));

for b = 1:numel(bc_vec)
  a_lim_vec_vanilla(b) = ALimMulitplexFun(S0, B0, RNP1, RNP2, bc_vec(b), k, wc,rcg,1,kc) /mol/ 1e-15 /reaction_vol;
end

vanilla = figure;
cmap1 = brewermap(10,'Blues');
cmap2 = brewermap(10,'Greens');
cmap3 = brewermap(10,'Blues');
hold on



ylim([1e-4 1e3])

% patch([bc_vec fliplr(bc_vec)], [a_lim_vec_vanilla min(ylim)*ones(size(bc_vec))], cmap3(9,:),'FaceAlpha',.2)        % Below Lower Curve
patch([bc_vec fliplr(bc_vec)], [a_lim_vec_vanilla max(ylim)*ones(size(bc_vec))], cmap2(5,:),'FaceAlpha',.2)  

plot(bc_vec,a_lim_vec_vanilla,'Color',cmap1(5,:),'LineWidth',2)

grid on
xlabel('relative Cas13 background')
ylabel('RNA detection limit (fM)')

set(gca,'YScale','log')
set(gca,'XScale','log')
set(gca,'Fontsize',12)

saveas(vanilla,[FigPath 'vanilla_ncr.png'])

%%
close all
% make figure showing impact of multiplexing
bc_vec = logspace(-10,0);
n_vec = [1 10 100 1000];
a_lim_array_n = NaN(numel(n_vec),numel(bc_vec));

for n = 1:numel(n_vec)
  for b = 1:numel(bc_vec)
    a_lim_array_n(n,b) = ALimMulitplexFun(S0, B0, 2*sqrt(n_vec(n))*RNP1, RNP2, bc_vec(b), k, wc,rcg,n_vec(n),kc) /mol/ 1e-15 /reaction_vol;
  end
end

multiplex = figure;
hold on
lgd_str = {};
ylim([1e-4 1e3])

% Below Lower Curve
% patch([bc_vec fliplr(bc_vec)], [a_lim_array_n(end,:) min(ylim)*ones(size(bc_vec))], cmap3(7,:),'FaceAlpha',.2) 

for n = 1:numel(n_vec)
  patch([bc_vec fliplr(bc_vec)], [a_lim_array_n(n,:) max(ylim)*ones(size(bc_vec))], cmap2(5,:),'FaceAlpha',.15)  
end
p = [];
for n = 1:numel(n_vec)
  p = [p plot(bc_vec,a_lim_array_n(n,:),'Color',cmap1(5+n,:),'LineWidth',2)];
  lgd_str = [lgd_str{:} {['NCR (' num2str(n_vec(n)) 'x)']}];
end

grid on
xlabel('background Cas13  nuclease actvity')
ylabel('RNA detection limit (fM)')
legend(p,lgd_str{:},'Location','southeast')
set(gca,'YScale','log')
set(gca,'XScale','log')
set(gca,'Fontsize',12)

saveas(multiplex,[FigPath 'multiplex_ncr.png'])


%% make figure showing impact of of RNP1
close all
bc = 1e-6;
rnp1_vec = logspace(-2,4)*nanomol*reaction_vol; 
a_lim_array_rnp1 = NaN(numel(rnp1_vec),numel(n_vec));
a_lim_array_rnp1_multi = NaN(numel(rnp1_vec),numel(n_vec));

for r = 1:numel(rnp1_vec)
  for n = 1:numel(n_vec)
    a_lim_array_rnp1(r,n) = ALimMulitplexFun(S0, B0, rnp1_vec(r), RNP2, bc, k, wc,rcg,n_vec(n),kc) /mol/ 1e-15 /reaction_vol;    
  end
end
rnp1_plot = rnp1_vec/nanomol/reaction_vol;

rnp1_v1 = figure;
cmap2 = brewermap(10,'Greens');
hold on
ylim([10^-2.5 1e3])
% patch([rnp1_plot fliplr(rnp1_plot)], [a_lim_array_rnp1(:,end)' min(ylim)*ones(size(bc_vec))], cmap3(7,:),'FaceAlpha',.2) 
for n = 1:numel(n_vec)
  patch([rnp1_plot fliplr(rnp1_plot)], [a_lim_array_rnp1(:,n)' max(ylim)*ones(size(bc_vec))], cmap2(5,:),'FaceAlpha',.15)  
end

p = [];
for n = 1:numel(n_vec)
  p = [p plot(rnp1_plot,a_lim_array_rnp1(:,n),'Color',cmap1(5+n,:),'LineWidth',2)];
end
grid on
xlabel('RNP1 concentration (nM)')
ylabel('RNA detection limit (fM)')
legend(p,lgd_str{:},'Location','southwest')
set(gca,'YScale','log')
set(gca,'XScale','log')
set(gca,'Fontsize',12)
% ylim([1e-4 1e3])
saveas(rnp1_v1,[FigPath 'rnp1_ncr.png'])