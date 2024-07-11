function y = equal2(a,b)
%EQUAL2 此处显示有关此函数的摘要
%   此处显示详细说明
if max(abs(a-b))<1e-10
    y=true;
else
    y=false;
end
end

