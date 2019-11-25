%% Artifial network
clear;
clc;
dir = 'D:\Ximing\Study\RESEARCH\UPenn\Source_Networks\';
%File = '100466178325794757407.edges';
%File = '117798157258572080176.edges';
%File = '113455290791279442483.edges';
%File = '110581012109008817546.edges';
File = '103241736833663734962.edges';
%File = 'yeastinter_st.txt';
%File = 'prisoninter_st.txt';
%File = 'opsahl-usairport.txt';
%File = 'moreno_sheep_sheep.txt';
%File = 'moreno_rhesus_rhesus.txt';
%File = 'moreno_oz_oz.txt';
%File = 'moreno_mac_mac.txt';
%File = 'moreno_highschool_highschool.txt';
%File = 'moreno_health_health.txt';
%File = 'leader2inter_st.txt';
File = 'Wiki-Vote.txt';
%File = 'foodweb-baydry.txt';
%File = 'coli1_1Inter_st.txt';
%A = RdSocialNetwork([dir File]);
A = RdNetwork([dir File]);
A(A~=0) = 1;    % remove the weight
A = A - diag(diag(A));  % remove diagonal

%% Random network
n = 500;
p = log(n)/n;
A = rand(n,n) < p;
A = A - diag(diag(A));

%% generate random chung-lu graph
n = 1500;
isDigraph = 1;
isSimple = 1;

beta = 5;
Delta = 120;
d = 40;
w = PowerLawCoef(beta, Delta, d, n);
w_cl = [w, w];
A = GenerateGraphs(n, w_cl, 'CL', isDigraph, isSimple);
%%
hh = hist(sum(A), 30);
set(gca, 'fontsize', 30);
hx = xlabel('In-degrees of vertices in the digraph $G$');
set(hx,'Interpreter','latex','fontsize',40);
hy = ylabel('Number of vertices');
set(hy,'Interpreter','latex','fontsize',40);
%% Display information
rho_A = max(abs(eig(A)));
rho_sym_A = max(eig( A+A'))/2;
w_A = max(abs(imag(eig(A))));
B = (A - A')/(2*1i);
w_B = max(abs(eig(B)));
disp(['Reciprocity of A: ', num2str(trace(A^2)/sum(sum(A)))]); 

%% Moment estimation framework
for r = 2:10
    shape = 'circle';
    optimize = 1;
    disp(['order #: ', num2str(r)]);
    Info(r-1).circle = MomentEstimationFramework(A, r, shape, optimize, 0);
    shape = 'square';
    Info(r-1).square = MomentEstimationFramework(A, r, shape, optimize, 0);
end
%% Single estimation
shape = 'circle';
optimize = 1;
r = 6;
disp(['order #: ', num2str(r)]);
out_circle = MomentEstimationFramework(A, r, shape, optimize, 0);
shape = 'square';
out_square = MomentEstimationFramework(A, r, shape, optimize, 0);


%% symmetrized upper bound
ratio_sym = [];
for order = 2:10
    ratio_sym = [ratio_sym; computeRhoSym(A, order)/rho_A];
end
%%
ratio_upper = [];
ratio_lower = [];
ratio_upper_refine = [];
ratio_lower_refine = [];
for i = 1:9
    rho_low = max(Info(i).circle.rho_low, Info(i).square.rho_low);
    rho_low_refined = max([rho_low, Info(i).circle.rho_low_refined, Info(i).square.rho_low_refined]);
    ratio_lower = [ratio_lower; rho_low/rho_A];
    ratio_lower_refine = [ratio_lower_refine; rho_low_refined/rho_A];
    rho_upp = min(Info(i).circle.rho_upp, Info(i).square.rho_upp);
    rho_upp_refined = min([rho_upp, Info(i).circle.rho_upp_refined, Info(i).square.rho_upp_refined]);
    ratio_upper = [ratio_upper; rho_upp/rho_A];
    ratio_upper_refine = [ratio_upper_refine; rho_upp_refined/rho_A];
end
%%
figure();
hold on;
plot(2:10, ratio_lower, 'bs','markerfacecolor','b','markersize',15);
hlower = plot(2:10, ratio_lower_refine,'b','linewidth',5);
hupper = plot(2:10, ratio_upper, 'r','linewidth',5);
plot(2:10, ratio_upper, 'rs','markerfacecolor','r','markersize', 15);

hupperrefine = plot(2:10, ratio_upper_refine,'g','linewidth',5);
plot(2:10, ratio_upper_refine,'gs','markerfacecolor','g','markersize', 15);
%plot(2:10, ratio_sym,'y','linewidth',3);
%plot(2:10, ratio_sym,'ys','markerfacecolor','y','markersize', 10);
hold off;
hall = [hlower, hupper, hupperrefine];
leg1 = legend(hall,'Lower bound $\underline{\rho}_r^\star$','Upper bound $\overline{\rho}_r^\star$','Refined upper bound $\overline{\varrho}_r^\star$');
set(leg1,'Interpreter','latex');
set(leg1,'FontSize',30);
hx = xlabel('Order of subgraphs considered');
set(hx, 'Interpreter','latex','fontsize',40);
hy = ylabel('Normalized bounds on $\lambda_n$');
set(hy, 'Interpreter','latex','fontsize',40);
set(gca,'fontsize',30);
axis([1.9, 10.1, -0.01, 10.1]);

