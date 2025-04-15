# A collection of functions associated with immersed network data types. 
# This is where force calculations, spread and interp, as well as 
import HDF5



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

function MoveNetwork!(mymesh::IBMProject.LagMesh,x::Vector{T},y::Vector{T}) where T <: Real
	# First lets have some sanity checks to make sure that 
	if ~(size(x) == size(y))
		throw(ArgumentError("Size of displacements does not match"));
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


function ElasticForceCalc(mymesh::LagMesh,lame1::T,lame2::T) where T <: Real
	Force = [[0.0,0.0] for i in 1:mymesh.Mv]

	for i = 1:mymesh.Fv

	end
end

