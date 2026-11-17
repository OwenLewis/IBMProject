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
	for v in mymesh.vertices
		v.area = 0.0;
	end


	#Some Auxilliar quantities to check for inverted faces
	ninverted = 0;
	firstbad = 0;
	worstJ = Inf;

	#First we loop through the faces
	for i in 1:mymesh.Mf
		#Who are this face's vertices?
		vert1,vert2,vert3 = mymesh.faces[i].vertices;
		#Update this face's displacement tensor
		mymesh.faces[i].Xmat = TriangleEdgeMatrix(vert1.coords, vert2.coords, vert3.coords);
		#Update this face's center
		mymesh.faces[i].centroid = TriangleCentroid(vert1.coords,vert2.coords,vert3.coords);
		#Calculate this face's area
		area = TriangleArea(vert1.coords,vert2.coords,vert3.coords);
		mymesh.faces[i].area = area;
		# And give 1/3 to every one of the vertices
		vert1.area += area/3;
		vert2.area += area/3;
		vert3.area += area/3;

		#Now lets check to see if the face has been inverted. The area
		#Calculation would have squashed this with abs() call. 
		J = LinearAlgebra.det(mymesh.faces[i].Xmat) / LinearAlgebra.det(mymesh.faces[i].Smat)
		if J <= 0
			ninverted += 1;
			if firstbad == 0
				firstbad = i
			end
			worstJ = min(worstJ, J);
		end

	end

	#Throw warnings if any face was inverted
	if ninverted > 0
		@warn "Mesh has inverted triangles" count=ninverted of=mymesh.Mf first_face=firstbad min_J=worstJ maxlog=10
	end

	#Finally we loop through the edges
	for i in 1:mymesh.Me
		#Who are this edge's vertices?
		vert1,vert2 = mymesh.edges[i].vertices;
		#Calcualte this edge's center point and assign. 
		mymesh.edges[i].center = Midpoint(vert1.coords,vert2.coords);
	end
end


# This function moves a network by a given velocity field. 
# It comes in three types
# Displacements as separate arrays
function MoveNetwork!(mymesh::LagMesh,x::Vector{<: Real},y::Vector{<: Real})
	# First lets have some sanity checks to make sure that 
	if ~(size(x) == size(y))
		throw(ArgumentError("Size of displacements do not match"));
	end
	if ~(length(x) == mymesh.Mv)
		throw(ArgumentError("Size of displacements does not match mesh vertices"));
	end

	# Update the coordinates of all the veritices 
	for i = 1:mymesh.Mv
		mymesh.vertices[i].coords += Point2(x[i],y[i]);
	end
	# Then update the mesh
	DefUpdate!(mymesh);
end

# Displacements as a vector of vectors
function MoveNetwork!(mymesh::LagMesh,X::Vector{<: Vector{<: Real}})
	# First lets have some sanity checks to make sure that 
	if ~(length(X[1]) == 2)
		throw(ArgumentError("Size of displacements is not appropriate"));
	end
	if ~(length(X) == mymesh.Mv)
		throw(ArgumentError("Size of displacements does not match mesh vertices"));
	end

	# Update the coordinates of all the veritices 
	for i = 1:mymesh.Mv
		mymesh.vertices[i].coords += Point2(X[i]);
	end
	# Then update the mesh
	DefUpdate!(mymesh);
end

# Displacements as a vector of 2D vectors
function MoveNetwork!(mymesh::LagMesh,X::Vector{<:StaticArrays.SVector{2,<:Real}})
	# First lets have some sanity checks to make sure that 
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


function ElasticForceCalc(mymesh::LagMesh,mu::T1,kappa::T2) where {T1 <: Real,T2 <: Real}
	#This calculates the elastic force in the LagMesh due to its current configuration
	# mu is the shear modulus 
	# kappa is the bulk modulus (elastic resistance to area increases)
	# Allocate an array of 2D points, currently everyone is set to zero.
	ForceDensity = zeros(Point2, mymesh.Mv)
	# A = StaticArrays.SMatrix{2,2,Float64}(0,0,0,0)

	for i = 1:mymesh.Mf
		
		#Grab reference to the face we're on. 
		face = mymesh.faces[i]

		#Now we calculate the deformation gradient tensor
		S = face.Smat;
		A = face.Xmat*face.Sinv;
		#And its associated quantities
		Ainv = LinearAlgebra.inv(A);
		AT = LinearAlgebra.transpose(A);
		AinvT = LinearAlgebra.transpose(Ainv);
		J = LinearAlgebra.det(A);

		#We also need the determinant of the original matrix
		Delta = LinearAlgebra.det(S);

		EnergyGradient = mu.*(A./J - LinearAlgebra.tr(A*AT).*AinvT./(2*J)) + kappa.*(J.- 1.0).*J.*AinvT;


		#Now we loop through the 3 points. Each has its own dA/dX
		#The calculations are what SHOULD appear in Strychalski, Copos, et al. 
		#One day I need to double check these things. 
		dAdX = (StaticArrays.SVector{2,Float64}((S[2,1]-S[2,2])/Delta,(S[1,2]-S[1,1])/Delta),
				StaticArrays.SVector{2,Float64}(S[2,2]/Delta,-S[1,2]/Delta),
				StaticArrays.SVector{2,Float64}(-S[2,1]/Delta,S[1,1]/Delta));

		#Loop through the vertices around this face. 
		#Multiply force density by face area, and then divide by vertex area
		for (k,vertex) in enumerate(face.vertices)
 			Force = -(EnergyGradient*dAdX[k]).*face.ref_area;
 			ForceDensity[vertex.id] += Force./vertex.ref_area;


		end


	end


	return ForceDensity
end

