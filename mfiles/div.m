function divout = div(space,var,direction,divType)
%This fucntion evaluates the derivative at j faces
if strcmp(direction,'i') == 1 && strcmp(divType,'central') == 1
    divout = icenterdiv(space,var);
elseif strcmp(direction,'j') == 1 && strcmp(divType,'central') == 1
    divout = jcenterdiv(space,var);
elseif strcmp(direction,'k') == 1 && strcmp(divType,'central') == 1
    divout = kcenterdiv(space,var);
elseif strcmp(direction,'i') == 1 && strcmp(divType,'old') == 1
    divout = iolddiv(space,var);
elseif strcmp(direction,'k') == 1 && strcmp(divType,'old') == 1
    divout = kolddiv(space,var);
elseif strcmp(direction,'i') == 1 && strcmp(divType,'fft') == 1
    divout = ifftdiv(space,var);
elseif strcmp(direction,'k') == 1 && strcmp(divType,'fft') == 1
    divout = kfftdiv(space,var);
elseif strcmp(direction,'i') == 1 && strcmp(divType,'5point') == 1
    divout = ifivepoint(space,var);
elseif strcmp(direction,'j') == 1 && strcmp(divType,'5point') == 1
    divout = jfivepoint(space,var);
elseif strcmp(direction,'k') == 1 && strcmp(divType,'5point') == 1
    divout = kfivepoint(space,var);
elseif strcmp(direction,'i') == 1 && strcmp(divType,'7point') == 1
    divout = isevenpoint(space,var);
elseif strcmp(direction,'j') == 1 && strcmp(divType,'7point') == 1
    divout = jsevenpoint(space,var);
elseif strcmp(direction,'k') == 1 && strcmp(divType,'7point') == 1
    divout = ksevenpoint(space,var);
end

end

function out = icenterdiv(space,var)
[m,n,p] = size(space);
out = zeros(size(var));
for i = 2:m-1
    for j= 2:n-1
        for k = 2:p-1
            out(i,j,k) = (var(i+1,j,k)+var(i+1,j-1,k)-var(i-1,j,k)-var(i-1,j-1,k))/...
                (space(i+1,j,k)+space(i+1,j-1,k)-space(i-1,j,k)-space(i-1,j-1,k));
        end
    end
end

end

function out = kcenterdiv(space,var)
[m,n,p] = size(space);
out = zeros(size(var));
for i = 2:m-1
    for j= 2:n-1
        for k = 2:p-1
            out(i,j,k) = (var(i,j,k+1)+var(i,j-1,k+1)-var(i,j,k-1)-var(i,j-1,k-1))/...
                (space(i,j,k+1)+space(i,j-1,k+1)-space(i,j,k-1)-space(i,j-1,k-1));
        end
    end
end

end

function out = jcenterdiv(space,var)
[m,n,p] = size(space);
out = zeros(size(var));
for i = 2:m-1
    for j= 2:n-1
        for k = 2:p-1
            out(i,j,k) = (var(i,j,k)-var(i,j-1,k))/...
                (space(i,j,k)-space(i,j-1,k));
        end
    end
end
end

function out = iolddiv(space,var)
[m,n,p] = size(space);
out = zeros(size(var));
for i = 2:m-1
    for j= 2:n-1
        for k = 2:p-1
            out(i,j,k) = (var(i,j,k)-var(i-1,j,k))/...
                (space(i,j,k)-space(i-1,j,k));
        end
    end
end
end

function out = kolddiv(space,var)
[m,n,p] = size(space);
out = zeros(size(var));
for i = 2:m-1
    for j= 2:n-1
        for k = 2:p-1
            out(i,j,k) = (var(i,j,k)-var(i,j,k-1))/...
                (space(i,j,k)-space(i,j,k-1));
        end
    end
end
end

function out = ifftdiv(space,var)
[~,n,p] = size(space);
L = space(end,1,1)-space(1,1,1);
out = zeros(size(var));
for j = 1:n
    for k = 1:p
        out(:,j,k) = fftdiff(var(:,j,k),L,1);
    end
end
end

function out = kfftdiv(space,var)
[m,n,~] = size(space);
W = space(1,1,end)-space(1,1,1);
out = zeros(size(var));
for i = 1:m
    for j = 1:n
        out(i,j,:) = fftdiff(squeeze(var(i,j,:)),W,1);
    end
end
end

function out = ifivepoint(space,var)
[m,n,p] = size(space);
out = zeros(size(var));

c = zeros(m,5);
for i = 3:m-2
    c(i,:) = fdcoeffF(1,space(i,1,1),space(i-2:i+2,1,1));
end

for i = 3:m-2
    for j = 1:n
        for k = 1:p
            out(i,j,k) = c(i,:)*var(i-2:i+2,j,k);
        end
    end
end
end

function out = jfivepoint(space,var)
[m,n,p] = size(space);
out = zeros(size(var));

c = zeros(n,5);
for j = 3:n-2
    c(j,:) = fdcoeffF(1,space(1,j,1),space(1,j-2:j+2,1));
end

for i = 1:m
    for j = 3:n-2
        for k = 1:p
            out(i,j,k) = c(j,:)*var(i,j-2:j+2,k)';
        end
    end
end
end

function out = kfivepoint(space,var)
[m,n,p] = size(space);
out = zeros(size(var));

c = zeros(p,5);
for k = 3:p-2
    c(k,:) = fdcoeffF(1,space(1,1,k),space(1,1,k-2:k+2));
end

for i = 1:m
    for j = 1:n
        for k = 3:p-2
            out(i,j,k) = c(k,:)*squeeze(var(i,j,k-2:k+2));
        end
    end
end
end

function out = jsevenpoint(space,var)
[m,n,p] = size(space);
out = zeros(size(var));
for i = 1:m
    for j = 5:n-4
        for k = 1:p
            c = fdcoeffF(1,space(i,j,k),space(i,j-3:j+3,k));
            out(i,j,k) = c*var(i,j-3:j+3,k)';
        end
    end
end
end
