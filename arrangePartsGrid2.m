function [ parts ] = arrangePartsGrid2( N,L )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
spacing = L/N;
x = -L/2+spacing/2;
y = -L/2+spacing/2;
parts = zeros(2*N*N,1);
%parts = zeros(N,2);
for i = 1:N
    for j = 1:N
         parts(2*((i-1)*N+j)-1) = x;
         parts(2*((i-1)*N+j)) = y;
%        parts(i,1) = x;
%         parts(i,2) = y;
        y = y + spacing;
    end
    y = -L/2+spacing/2;
    x = x+spacing;
end
        

end

