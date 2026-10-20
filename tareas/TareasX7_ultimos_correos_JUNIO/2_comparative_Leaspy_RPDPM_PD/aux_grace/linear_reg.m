function yCalc=linear_reg(x,y)
X=[ones(length(x),1),x];
b=X\y;
yCalc = X*b; 
end