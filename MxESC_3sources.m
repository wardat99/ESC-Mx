clear
close all
clc
warning off;
addpath('network','functions');

load('myThreeSources.mat')

mu = 1;                    %best at   1 nmi=0.7468 
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