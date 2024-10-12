
function [mu,sd,med,ent,skew,kurt]=texturefeature(yyy)
%%%%%%%%% TEXTURE ANALYSIS %%%%%%%%%
%%%1.MEAN%%%%%
mu=mean2(yyy);

%%%2.STD %%%%%
sd=std(yyy);

%%%3.MEDIAN%%%%%
med=median(yyy);

%%%4.ENTROPY %%%%%
ent=entropy(yyy);

%%%5.SKEWNESS %%%%%
skew=skewness(yyy);

%%%6.KURTOSIS %%%%%
kurt=kurtosis(yyy);
end