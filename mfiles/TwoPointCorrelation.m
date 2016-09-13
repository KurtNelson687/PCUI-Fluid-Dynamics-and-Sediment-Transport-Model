function [TwoPoint] = TwoPointCorrelation(var,r_all,heightIndex,plane,direction)
% Kurt NElson; 4/8/2-16
% 2PointCorrelation - Returns the 2PointCorrelation of var for all points in the
% domain seperated by distance(s) r. The plane the correlation is taken
% over is given by plane ( in PCUI 13 = x-z plave, 2 = vertical, 3 = cross
% stream). heightIndex gives the index inwhich things are summed across.
[m,n,p] = size(var);
numerator = 0;
denominator = 0;
TwoPoint = zeros(size(r_all));

if plane == 13 && direction == 1; %aveage
    for r = 1:length(r_all)
        for i = 1:m
            for k = 1:p
                rtemp =r+i;
                if r+i>m
                    rtemp = r+i-m;
                end
                numerator = numerator+var(i,heightIndex,k)*var(rtemp,heightIndex,k);
                denominator = denominator+var(i,heightIndex,k)^2;
            end
        end
        numerator = numerator/(m+p);
        denominator = denominator/(m+p);
        TwoPoint(r)= numerator/denominator;
        numerator= 0;
        denominator = 0;
    end
elseif plane == 13 && direction == 3;
    for r = 1:length(r_all)
        for i = 1:m
            for k = 1:p
                rtemp =r+k;
                if r+k>p
                    rtemp = r+k-p;
                end
                numerator = numerator+var(i,heightIndex,k)*var(i,heightIndex,rtemp);
                denominator = denominator+var(i,heightIndex,k)^2;
            end
        end
        numerator = numerator/(m+p);
        denominator = denominator/(m+p);
        TwoPoint(r)= numerator/denominator;
        numerator= 0;
        denominator = 0;
    end
end

end

