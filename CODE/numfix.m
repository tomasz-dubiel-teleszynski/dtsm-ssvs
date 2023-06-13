function iw = numfix(iw)

if max(iw) < -700
    const = 350 - max(iw);
    iw = iw + const;
end
if max(iw) > 350
    const = max(iw) - 350;
    iw = iw - const;
end
                
end