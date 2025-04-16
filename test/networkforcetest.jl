
using IBMProject
using GLMakie

flag = 2

if flag == 0
	mymesh = ReadNetwork("smallgrid.h5")
elseif flag == 1
	mymesh = ReadNetwork("medgrid.h5");
elseif flag == 2
	mymesh = ReadNetwork("largegrid.h5")
else
	points = [[0.0,0.0],[1.0,0.0],[0.0,1.0]]
	faces = [[1,2,3]]
	mymesh = LagMesh(points,faces);
end







println("Grid Constructed")

println("Final grid has:")
println("	",mymesh.Mv," vertices,")
println("	",mymesh.Me," edges,")
println("	",mymesh.Mf," faces.")


U = [x.coords[1]*0.2 for x in mymesh.vertices];
V = [x.coords[2]*0.2 for x in mymesh.vertices];


println("Moving the grid now")
IBMProject.MoveNetwork!(mymesh,U,V);
X = [x.coords[1] for x in mymesh.vertices];
Y = [x.coords[2] for x in mymesh.vertices];


println("Calculating force")
forcedensity = IBMProject.ElasticForceCalc(mymesh,1,0)
hcomp = [x[1] for x in forcedensity];
vcomp = [x[2] for x in forcedensity];

strength = [sqrt.(x[1].^2 + x[2].^2) for x in forcedensity];

fig = Figure()
ax = Axis(fig[1,1])
foo = arrows!(ax,X,Y,hcomp,vcomp, lengthscale = 0.03, arrowcolor = strength, linecolor = strength,linewidth = 3)
for i = 1:length(mymesh.edges)
	xs = [x.coords[1] for x in mymesh.edges[i].vertices];
	ys = [x.coords[2] for x in mymesh.edges[i].vertices];
	lines!(ax,xs,ys,color = :black,linewidth = 0.7)
end

fig
