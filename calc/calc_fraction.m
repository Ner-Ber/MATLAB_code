function [Numerator,Denominator]=calc_fraction(n,maxDenominator)

Denominator=2:maxDenominator;
Numerator =1:floor(maxDenominator/2);

tmp=Numerator./Denominator;

[~,index]=min(abs(tmp-n));

Numerator=Numerator(index);
Denominator=Denominator(index);
