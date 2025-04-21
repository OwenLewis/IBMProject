# A collection of functions associated with immersed network data types. 
# This is where force calculations, spread and interp, as well as 
import HDF5
import StaticArrays
import LinearAlgebra


# This is a simple function to read a network in from an HDF5 file
function ReadNetwork(filename::String)
	fid = HDF5.h5open(filename)

	readpoints = HDF5.read(fid,"p");
	points = [readpoints[i,:] for i in 1:size(readpoints,1)];

	readfaces = Int.(read(fid,"t"));
	faces = [readfaces[i,:] for i in 1:size(readfaces,1)];
	mymesh = LagMesh(points,faces);
	close(fid)

	return mymesh
end


# These functions update all of the eulerian characteristics of a network. 
# The assumption is that all of the point coords have been updated so that 
# They are no longer equal to the ref_coords
function DefUpdate!(mymesh::LagMesh)
	[x.area = 0.0 for x in mymesh.vertices];
	#First we loop through the faces
	for i in 1:mymesh.Mf
		#Who are this face's vertices?
		one,two,three = mymesh.faces[i].vertices;
		#Update this face's displacement tensor
		mymesh.faces[i].Xmat = hcat(two.coords .- one.coords,three.coords .- one.coords);
		#Update this face's center
		mymesh.faces[i].centroid = (one.coords + two.coords + three.coords)/3
		#Calculate this face's area
		area = abs(one.coords[1]*(two.coords[2] - three.coords[2]) + 
					two.coords[1]*(three.coords[2] - one.coords[2]) + 
					three.coords[1]*(one.coords[2] - two.coords[2]))/2;
		mymesh.faces[i].area = area;
		# And give 1/3 to every one of the vertices
		one.area += area/3;
		two.area += area/3;
		three.area += area/3;
	end

	#Finally we loop through the edges
	for i in 1:mymesh.Me
		#Who are this edge's vertices?
		one,two = mymesh.edges[i].vertices;
		#Calcualte this edge's center point and assign. 
		mymesh.edges[i].center = (one.coords + two.coords)./2
	end
end


# This function moves a network by a given velocity field. 
# It comes in two types

function MoveNetwork!(mymesh::LagMesh,x::Vector{T},y::Vector{T}) where T <: Real
	# First lets have some sanity checks to make sure that 
	if ~(size(x) == size(y))
		throw(ArgumentError("Size of displacements do not match"));
	end
	if ~(length(x) == mymesh.Mv)
		throw(ArgumentError("Size of displacements does not match mesh vertices"));
	end

	# Update the coordinates of all the veritices 
	for i = 1:mymesh.Mv
		mymesh.vertices[i].coords += [x[i],y[i]];
	end
	# Then update the mesh
	DefUpdate!(mymesh);
end


function MoveNetwork!(mymesh::LagMesh,X::Vector{Vector{T}}) where T <: Real
	# First lets have some sanity checks to make sure that 
	if ~(length(X[1]) == 2)
		throw(ArgumentError("Size of displacements is not appropriate"));
	end
	if ~(length(X) == mymesh.Mv)
		throw(ArgumentError("Size of displacements does not match mesh vertices"));
	end

	# Update the coordinates of all the veritices 
	for i = 1:mymesh.Mv
		mymesh.vertices[i].coords += X[i];
	end
	# Then update the mesh
	DefUpdate!(mymesh);
end


function ElasticForceCalc(mymesh::LagMesh,mu::T,kappa::T) where T <: Real
	#This calculates the elastic force in the LagMesh due to its current configuration
	# mu is the shear modulus 
	# kappa is the bulk modulus (elastic resistance to area increases)
	ForceDensity = [[0.0,0.0] for i in 1:mymesh.Mv]
	A = StaticArrays.SMatrix{2,2,Float64}(0,0,0,0)

	for i = 1:mymesh.Mf
		# I'm indexing the points zero, one, two to correspond with Strychalski, Copos, et al
		# I'm only doing that here. 
		zero,one,two = mymesh.faces[i].vertices;
		ind0 = findfirst(x -> x==zero,mymesh.vertices);
		ind1 = findfirst(x -> x==one,mymesh.vertices);
		ind2 = findfirst(x -> x==two,mymesh.vertices);

		#Now we calculate the deformation gradient tensor
		S = mymesh.faces[i].Smat;
		A = mymesh.faces[i].Xmat*mymesh.faces[i].Sinv;
		#And its associated quantities
		Ainv = LinearAlgebra.inv(A);
		AT = LinearAlgebra.transpose(A);
		AinvT = LinearAlgebra.transpose(Ainv);
		J = LinearAlgebra.det(A);

		#We also need the determinant of the original matrix
		Delta = LinearAlgebra.det(S);

		EnergyGradient = mu.*(A./J - LinearAlgebra.tr(A*AT).*AinvT./(2*J)) + kappa.*(J.- 1.0).*J.*AinvT;


		#Now we loop through the 3 points. Each has its own dA/dX
		#Point zero
		dAdX = StaticArrays.SVector{2,Float64}((S[2,1]-S[2,2])/Delta,(S[1,2]-S[1,1])/Delta)
		Force = -(EnergyGradient*dAdX).*mymesh.faces[i].ref_area;
		ForceDensity[ind0] += Force./zero.ref_area
		#Point one
		dAdX = StaticArrays.SVector{2,Float64}(S[2,2]/Delta,-S[1,2]/Delta)
		Force = -(EnergyGradient*dAdX).*mymesh.faces[i].ref_area;
		ForceDensity[ind1] += Force./one.ref_area
		#Point two
		dAdX = StaticArrays.SVector{2,Float64}(-S[2,1]/Delta,S[1,1]/Delta)
		Force = -(EnergyGradient*dAdX).*mymesh.faces[i].ref_area;
		ForceDensity[ind2] += Force./two.ref_area

	end


	return ForceDensity
end

