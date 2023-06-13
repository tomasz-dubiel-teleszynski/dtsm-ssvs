function fi = bsESSnrat(fi1,logiw,logwnr,crit)

NMAX = 1000;
a = fi1;
b = 1;
N = 1;
while N <= NMAX
    c = (a + b)/2;
    if ((b - a)/2 < eps) || (getESSnrat(c,fi1,logiw,logwnr) >= crit) 
        break;
    end
    if sign(getESSnrat(c,fi1,logiw,logwnr) - crit) == sign(getESSnrat(a,fi1,logiw,logwnr) - crit)
        a = c;
    else
        b = c;
    end
    N = N + 1;
end
fi = c;

end