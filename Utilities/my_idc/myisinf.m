function ind=isinf(x);

ind=~(finite(x) | isnan(x) | isstr(x));
