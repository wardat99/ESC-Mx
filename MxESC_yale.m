clear
close all
clc
warning off;
addpath('network','functions');

load('myYale.mat')

mu = 10;                   %best at 10 nmi=0.9830
lambda = 0.1;              %best at 0.1
pr=0.6;                    %best at 0.6   

for i=1:10
    fprintf('--------MxESC start, attempt number %d--------\n', i);
[Plabel,Timecost(i)] = MxESC(A,mu,lambda,numClust,pr);
   fprintf('--------MxESC end, attempt number %d--------\n', i);
acc(i) =  Compute_accuracy(truth,Plabel);
        [~, nmi(i), ~] = compute_nmi(truth,Plabel);
        [f(i),p(i),r(i)] = compute_f(truth,Plabel);
        if (min(truth)==0)
            [AR(i)]=RandIndex(truth+1,Plabel);
        else
            [AR(i)]=RandIndex(truth,Plabel);
        end 
end
fprintf('Acc: %.2f  (%.2f)\n' , mean(acc), std(acc));
fprintf('nmi: %.4f  (%.4f)\n' , mean(nmi), std(nmi));
fprintf('AR: %.4f   (%.4f)\n' , mean(AR), std(AR));
fprintf('F: %.4f    (%.4f)\n' , mean(f), std(f));
fprintf('P: %.4f    (%.4f)\n' , mean(p), std(p));
fprintf('R: %.4f    (%.4f)\n' , mean(r), std(r));
fprintf('Timecost: %.4f  \n\n' , mean(Timecost));

% mesh(nmi)
% surf(nmi)
% colormap(jet)
% xlabel('\gamma_2')
% ylabel('\gamma_1')
% zlabel('NMI value')
% xticklabels({'0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1'})
% yticklabels({'0','1','2','3','4','5','6','7','8','9','10'})