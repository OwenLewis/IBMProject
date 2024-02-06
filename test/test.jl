using IBMProject
using Plots

plotlyjs()

M = 500;
Nx = 500;
Ny = 500;

mygrid = PeriodicEulGrid(Nx,Ny);

theta = 0:2*pi/M:2*(pi-1/M);
x = @. 0.25*cos(theta) + 0.5;
y = @. 0.25*sin(theta) + 0.5;

mybnd = PeriodicLagBnd(x,y);

bnddata = ones(M);

griddata = IBMProject.ScalarBndSpread(bnddata,mybnd,mygrid);
A = heatmap(mygrid.X[:,1],mygrid.Y[1,:],griddata)


println("Test case:")
println("Integral on Boundary: ", IBMProject.BndIntegral(bnddata,mybnd))
println("Integral on Grid:     ", IBMProject.GridIntegral(griddata,mygrid))


as_struct = ScalarGridData(griddata,mygrid);

L = IBMProject.MakePeriodicLaplacian(mygrid); #Lets create a laplacian

#Now (by hand), we'll create an inverse laplacian.
foo = L.Eigenvalues;
bar = 1./foo;
bar[1,1] = 0;
Linv = IBMProject.SimplePeriodicDifferentialOperator(mygrid,bar);


foo = IBMProject.ApplySimpleOperator(as_struct,L);
bar = IBMProject.ApplySimpleOperator(as_struct,Linv);


B = heatmap(mygrid.X[:,1],mygrid.Y[1,:],foo.U);
C = heatmap(mygrid.X[:,1],mygrid.Y[1,:],bar.U);

A
