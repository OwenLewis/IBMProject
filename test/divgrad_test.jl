using IBMProject


L = 3;
H = 2;

testu(x,y) = @. x^2*(x-L)^2*y^2*(y-H)^2 - (H^5*L^5)/(900.0);
testux(x,y) = @. 2*x*(x-L)^2*y^2*(y-H)^2+ 2*x^2*(x-L)*y^2*(y-H)^2;
testuy(x,y) = @. 2*x^2*(x-L)^2*y*(y-H)^2 + 2*x^2*(x-L)^2*y^2*(y-H);

Nx = 1000;
Ny = 1000;

mygrid = PeriodicEulGrid(L,H,Nx,Ny);

u = testu(mygrid.X,mygrid.Y);
u_struct = ScalarGridData(u,mygrid);

gradux = testux(mygrid.X,mygrid.Y);
graduy = testuy(mygrid.X,mygrid.Y);

gradu_struct = VectorGridData(gradux,graduy,mygrid);
# nabla = MakePeriodicGradient(mygrid);
nabladot = IBMProject.MakePeriodicDivergence(mygrid);
nabla = IBMProject.MakePeriodicGradient(mygrid);

foou,fooy = ApplyGradientOperator(u,nabla);

foo_struct = ApplyGradientOperator(u_struct,nabla);

erroru = gradux .- foou;
errorv = graduy .- fooy;

error = @. sqrt(erroru^2 + errorv^2);

L1 = sum(error)*mygrid.dx*mygrid.dy
L2 = sqrt(sum(error.^2)*mygrid.dx*mygrid.dy)
L∞ = maximum(error)


println("L1 error of gradient: ",L1)
println("L2 error of gradient: ",L2)
println("Linf error of gradient: ",L∞)

U = u;
V = -U;

truediv = gradux - graduy;

vec_struct = VectorGridData(U,V,mygrid);

divstruct = ApplyDivergenceOperator(vec_struct,nabladot)

div = ApplyDivergenceOperator(U,V,nabladot);

error = abs.(truediv - div);

L1 = sum(error)*mygrid.dx*mygrid.dy;
L2 = sqrt(sum(error.^2)*mygrid.dx*mygrid.dy);
L∞ = maximum(error);

println("L1 error of gradient: ",L1)
println("L2 error of gradient: ",L2)
println("Linf error of gradient: ",L∞)