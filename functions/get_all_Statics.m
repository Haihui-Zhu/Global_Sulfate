
function [ta1,slope,b,R2,NRMSD,nmb,obsm,N] = get_all_Statics(X,Y)
nonnan = ~isnan(X) & ~isnan(Y);
X = X(nonnan);
Y = Y(nonnan);

if length(Y)<3
    fprintf('No enough input data for scatter\n')
    tN = sprintf('N = %d',numel(X));
    ta1=tN;
    slope=NaN;
    b=NaN;
    R2=NaN;
    nmb = NaN;
    N=NaN;
    obsm=NaN;
else
        % calc statistics
        R2 = corrcoef(X,Y,'rows','complete');
        R2 = R2(1,2)^2; % R2
        [slope, b] = organic_regress(Y,X);
        pm = '+-';
        obsm = mean(X);
        NRMSD = sqrt(mean((Y-X).^2))/obsm;
        nmb = mean(Y-X)/obsm;
        N = numel(X);

        tslope = sprintf('y = %.2fx%s%.1f',slope,pm((b<0)+1),abs(b));
        tRMSE = sprintf('NRMSD = %.2f',NRMSD);
        tr = ['R^2 = ' sprintf('%.2f',R2)];
        tN = sprintf('N = %d',N);

        ta1 = sprintf('%s\n%s\n%s\n%s',tr,tslope,tRMSE,tN);
end

end