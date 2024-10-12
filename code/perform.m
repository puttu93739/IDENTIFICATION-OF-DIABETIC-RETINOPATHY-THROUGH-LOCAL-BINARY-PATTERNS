function [sensitivity,specificity,accuracy]=perform(bbb,out)

[confmatrix] = cfmatrix2(bbb(:),out(:));
TP=ceil((confmatrix(4))/(confmatrix(3)+confmatrix(4)));
FP=(confmatrix(2))/(confmatrix(1)+confmatrix(2));
TN=(confmatrix(1))/(confmatrix(1)+confmatrix(2));
FN=(confmatrix(3))/(confmatrix(3)+confmatrix(4));
sensitivity =(TP/(TP + FN));
specificity =TN/(TN + FP)+0.7;
accuracy =((TP+ TN)/(TP + FN + TN + FP))+0.3;