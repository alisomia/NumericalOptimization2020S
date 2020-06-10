%This SCRIPT file stores some test functions for numerical experiment
%
% STRUCTURE
% -----------------------------------------------
% func [struct]
%       func.f : The Objective fuction.
%       func.g: The gradient function.
%       func.init : The initial point of the problem.
%
%----------------------------------------------
% Linting @ PKU
% 20200417

% Trigonometric Function

trig2.f = @(x) Trigonometric(x,1e2);
trig2.g = @(x) GradTrigonometric(x,1e2);
trig2.init = TrigonometricInit(1e2);
trig3.f = @(x) Trigonometric(x,1e3);
trig3.g = @(x) GradTrigonometric(x,1e3);
trig3.init = TrigonometricInit(1e3);
trig4.f = @(x) Trigonometric(x,1e4);
trig4.g = @(x) GradTrigonometric(x,1e4);
trig4.init = TrigonometricInit(1e4);


% tridiagonal function

trid2.f = @(x) Tridiagonal(x,1e2);
trid2.g = @(x) GradTridiagonal(x,1e2);
trid2.init =  TridiagonalInit(1e2);
trid3.f = @(x) Tridiagonal(x,1e3);
trid3.g = @(x) GradTridiagonal(x,1e3);
trid3.init =  TridiagonalInit(1e3);
trid4.f = @(x) Tridiagonal(x,1e4);
trid4.g = @(x) GradTridiagonal(x,1e4);
trid4.init =  TridiagonalInit(1e4);

% extended powell function
ep2.f = @(x) ExtendedPowell(x,1e2);
ep2.g = @(x) GradExtendedPowell(x,1e2);
ep2.init =  ExtendedPowellInit(1e2);
ep3.f = @(x) ExtendedPowell(x,1e3);
ep3.g = @(x) GradExtendedPowell(x,1e3);
ep3.init =  ExtendedPowellInit(1e3);
ep4.f = @(x) ExtendedPowell(x,1e4);
ep4.g = @(x) GradExtendedPowell(x,1e4);
ep4.init =  ExtendedPowellInit(1e4);

% mat sqrt problem
mat2.f = @(x) MatSqrt(x,10);
mat2.g = @(x) GradMatSqrt(x,10);
mat2.init = MatSqrtInit(10);
mat3.f = @(x) MatSqrt(x,32);
mat3.g = @(x) GradMatSqrt(x,32);
mat3.init = MatSqrtInit(32);
mat4.f = @(x) MatSqrt(x,100);
mat4.g = @(x) GradMatSqrt(x,100);
mat4.init = MatSqrtInit(100);



%%% SUB PROGRAM BEGIN
function fx = Trigonometric(x,n)
global cf
cf = cf+1;
 j = (1:n)';
%ffx = (i==j).*sin(x(j))+(i.*(i==j)+1).*cos(x(j));
%ffx = eye(n).*sin(x) + (i.*eye(n)+1).*cos(x);
%ffx = (sin(x))' + sum((i.*eye(n)+1).*cos(x),1);
ffx = j + n - ((sin(x))+ j.*cos(x) +sum(cos(x)));
%fx =sum((n+i-sum(ffx,1)).^2,2);
fx = sum(ffx.^2,1);
end



function gx = GradTrigonometric(x,n)
global cg
cg = cg+1;
%gx = zeros(size(x,1),1);
 j = (1:n)';
%ffx = (i==j).*sin(x(j))+(i.*(i==j)+1).*cos(x(j));
ffx = j+n-((sin(x))+ j.*cos(x) +sum(cos(x)));
%gfx =  diag(cos(x))-(diag(i)+1).*sin(x);
 %gx =-2*sum((n+i-sum(ffx,1)).*gfx,2);
 %gx = -2*sum((n+i-ffx).*gfx,2);
 gx = (2*j.*sin(x) -2*cos(x)).*ffx + 2*sum(ffx)*sin(x);
end
function x0 = TrigonometricInit(n)
    x0 = ones(n,1)./n;
end

function fx = ExtendedPowell(x,~)
global cf
cf = cf+1;
    x0 = x(4:4:end);
    x1 = x(3:4:end);
    x2 = x(2:4:end);
    x3 = x(1:4:end);
    y = (x3 + 10*x2).^2 + 5*(x1-x0).^2 + (x2-2*x1).^4 + 10*(x3 - x0).^4;
    fx = sum(y);
end

function gx = GradExtendedPowell(x,~)
global cg
cg = cg+1;
gx = zeros(size(x,1),1);
    x0 = x(4:4:end);
    x1 = x(3:4:end);
    x2 = x(2:4:end);
    x3 = x(1:4:end);
gx(4:4:end)  = -10*(x1-x0) - 40*(x3-x0).^3;
gx(3:4:end) = 10*(x1-x0) - 8*(x2-2*x1).^3;
gx(2:4:end) = 20*(x3+10*x2) + 4*(x2-2*x1).^3;
gx(1:4:end) = 2*(x3+10*x2) + 40*(x3-x0).^3;
end


function x = ExtendedPowellInit(n)
x = zeros(n,1);
x(4:4:end) = 3;
x(3:4:end) = 0;
x(2:4:end) = -1;
x(1:4:end) = 3;
end

function fx = Tridiagonal(x,n)
global cf
cf = cf+1;
    y = (2:n)'.*((2*x(2:n)-x(1:n-1)).^2);
    fx = sum(y);
end

function gx = GradTridiagonal(x,n)
global cg
cg = cg+1;
    gx  = zeros(size(x,1),1);
    gx(2:n) = 4*(2:n)'.*(2*x(2:n)-x(1:n-1));
    gx(1:n-1) = gx(1:n-1) - 2*((2:n)'.*(2*x(2:n)-x(1:n-1)));
end

function x = TridiagonalInit(n)
    x = ones(n,1);
end

function fx = MatSqrt(x,n)
global cf
cf = cf +1;
    y = sin((1:n^2));
    A = reshape(y,[n,n]);
    B =reshape(x,[n,n]);
    fx = norm(B*B-A*A,'fro')^2;
end

function gx = GradMatSqrt(x,n)
global cg
cg = cg +1;
    y = sin((1:n^2))';
    A = reshape(y,[n,n]);
    B =reshape(x,[n,n]);
    C = B*B-A*A;
    D = C*B'+B'*C;
    gx = 2*reshape(D,[n^2,1]);
end
function x = MatSqrtInit(n)
    x = 0.2*sin((1:n^2))';
end
%%% SUB PROGRAM END
    