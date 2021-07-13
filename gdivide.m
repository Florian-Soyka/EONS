function [answer] = gdivide(A,B)
%GADD Summary of this function goes here
%   Detailed explanation goes here
sizeA = size(A);
sizeB = size(B);
% diffRows=sizeA(1)-sizeB(1);
% diffColumns=sizeA(2)-sizeB(2);
% if(sum(sizeA) > 3 || sum(sizeB) > 3)
%     disp('Potential problem in gmultiply');
%     size(A)
%     size(B)
% end

answer = A./B;
end

