% for j = 1:3
%     for i = 1:5
% %         i
% %         j
%         "still here"
%         if i ==3
%             j
%             return
%         end
% %         if j ==2
% %             break
% %         end
%     end
% end
params.alpha = [1e-4,1e-3,1e-2,1,1e6]; %5
params.beta = [0.1,0.5,1,2]; %4
params.sigma_r = [0.01,0.5,1,2,5]; %5
params.gamma = [0.1,1,10,35]; %4
params.kappa = [1e-5,1e-4,1e-3,1e-2,1e-1,1,10,100,1e3,1e4,1e5]; %11
params.eta = [1e-1,1,10,1e2]; %4

alpha = params.alpha;
beta = params.beta;
sigma_r = params.sigma_r;
gamma = params.gamma;
kappa = params.kappa;
eta = params.eta;
count = 0;
for a = 2%length(alpha)
    for b = 3:length(beta)
        for s = 5:length(sigma_r)
            for g = 1:length(gamma)
                for k = 1:length(kappa)
                    for e = 1:length(eta)
                        [a,b,s,g,k,e]
                        count = count + 1;
                        if a == 3
                            return
                        end
                    end
                end
            end
        end
    end
end
