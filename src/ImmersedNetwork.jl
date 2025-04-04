abstract type AbstractVertex end
abstract type AbstractEdge end
abstract type AbstractFace end
abstract type AbstractNetwork end

abstract type AbstractVertexFunction end
abstract type AbstractEdgeFunction end
abstract type AbstractFaceFunction end
abstract type AbstractNetworkFunction end



#Lets define the struct for a mesh vertex
#And write it a couple of constructors


mutable struct MeshVertex <: AbstractVertex
	#Location of Vertex
	coords::Vector{Float64}
	ref_coords::Vector{Float64}

	#current and reference area
	area::Float64
	ref_area::Float64

	#Booleans that indicate if this is a boundary vertex and if it is linked to anything
	boundary::Bool
	linked::Bool
	#An integer to indicate boundary linkage (for later use)

	#Now we need vectors of this nodes faces and edges
	edges::Vector{AbstractEdge}
	faces::Vector{AbstractFace}


	#Elastic force (density?) at vertex (might want to move this to another type)
	# elastic_force::Vector{Float64}

end



function MeshVertex(X::Vector{T}) where T <: Real
	if ~(length(X) == 2)
		throw(ArgumentError("Location must be a 2D location"));
	end


	point::MeshVertex = MeshVertex(X,X,-1,-1,false,false,[],[]);
	return point

end


#Now lets define the struct for a mesh face
#And give it some construtors

mutable struct MeshFace <: AbstractFace
	#Location of face centroid
	centroid::Vector{Float64}
	ref_centroid::Vector{Float64}

	#current and reference area
	area::Float64
	ref_area::Float64

	#Booleans that indicate if this is a boundary face
	boundary::Bool

	#Now we need vectors of this nodes faces and vertices
	edges::Vector{AbstractEdge}
	vertices::Vector{AbstractVertex}

end

function MeshFace(one::MeshVertex,two::MeshVertex,three::MeshVertex)
	centroid = (one.coords + two.coords + three.coords)/3
	area = abs(one.coords[1]*(two.coords[2] - three.coords[2]) + 
				two.coords[1]*(three.coords[2] - one.coords[2]) + 
				three.coords[1]*(one.coords[2] - two.coords[2]))/2

	edge1 = MeshEdge(one,two)
	edge2 = MeshEdge(two,three)
	edge3 = MeshEdge(three,one)

	face::MeshFace = MeshFace(centroid,centroid,area,area,false,[edge1,edge2,edge3],[one,two,three])
	return face
end

function MeshFace(vertices::Vector{MeshVertex})
	if ~(length(vertices) == 3)
		throw(ArgumentError("Face must be defined by 3 vertices"));
	end
	one = vertices[1];
	two = vertices[2];
	three = vertices[3];
	face::MeshFace = MeshFace(one,two,three)
	return face
end


#Now lets define the struct for a mesh edge
#And give it some construtors

mutable struct MeshEdge <: AbstractEdge
	#Location of face centroid
	center::Vector{Float64}
	ref_center::Vector{Float64}

	#Booleans that indicate if this is a boundary face
	boundary::Bool

	#Now we need vectors of this nodes faces and vertices
	vertices::Vector{AbstractVertex}
	faces::Vector{AbstractFace}

end

function MeshEdge(one::MeshVertex,two::MeshVertex)
	center = (one.coords + two.coords)/2

	edge::MeshEdge = MeshEdge(center,center,false,[one,two],[])
	return edge
end



#Lets define an Lagrangian Mesh
mutable struct LagMesh <: AbstractNetwork
	#Number of vertices in the immersed network
	Mv::Int
	#Number of triangles in the immersed network
	Mf::Int

	#Location of immersed boundary points, and faces that connect them
	vertices::Vector{AbstractVertex}
	faces::Vector{AbstractFace}
	# edges::Vector{AbstractEdge}
end


function LagMesh(points::Vector{Vector{Float64}},triangles::Vector{Vector{S where S <: Int}})
	Mv = length(points)
	Mf = length(triangles)

	vertices = Vector{MeshVertex}(undef,Mv)
	for i = 1:Mv
		vertices[i] = MeshVertex(points[i])
	end

	faces = Vector{MeshFace}(undef,Mf)
	for i = 1:Mf
		ipoints = triangles[i]
		faces[i] = MeshFace(vertices[ipoints[1]],vertices[ipoints[2]],vertices[ipoints[3]])
	end

	mesh::LagMesh = LagMesh(Mv,Mf,vertices,faces)
	return mesh

end


function LinkMesh(mymesh::LagMesh)
	
end
# 	stretch_force::Vector{Vector{Float64}}
# 	# connect_force::Matrix{Float64} #This one won't be used until we have an immersed mesh. 
# 	#Arrays for constants that define constitutive laws for the above forces
# 	lame1::Vector{Float64}
# 	lame2::Vector{Float64}


# 	#Arrays which are used to flag and define points of the immersed network to which the bnd is connected
# 	# flag::Vector{Int}  #This one won't be used until we have an immersed mesh. 
# 	# connect_buddy{Int} #This one won't be used until we have an immersed mesh. 


# #We can make a periodic lagrangian boundary as long as we get two vectors of location points. 
# function PeriodicLagBnd(X::Vector{T},Y::Vector{T}) where T <: Real
# 	if ~(size(X) == size(Y))
# 		throw(ArgumentError("Size of inputs does not match"));
# 	end
# 	N = length(X);
# 	location = hcat(X,Y);
# 	taup = zeros(size(location));
# 	taum = zeros(size(location));
# 	norm0 = zeros(size(location));
	
# 	stretch_force = 0*location;

# 	kp = ones(N);
# 	km = ones(N);
# 	gammap = ones(N);
# 	gammam = ones(N);

# 	taup[1:N-1,:] = location[2:N,:] - location[1:N-1,:];
# 	taup[N,:] = location[1,:] - location[N,:];

# 	taum[2:N,:] = location[2:N,:] - location[1:N-1,:];
# 	taum[1,:] = location[1,:] - location[N,:];
	
# 	tau0 = (taup .+ taum)./2.0;


# 	dsp = @. sqrt(taup[:,1]^2 + taup[:,2]^2);
# 	dsm = @. sqrt(taum[:,1]^2 + taum[:,2]^2);
# 	ds0 = @. (dsp + dsm)/2;


# 	norm0[:,1] = -tau0[:,2];
# 	norm0[:,2] = tau0[:,1];

# 	bnd::PeriodicLagBnd = PeriodicLagBnd(N,location,stretch_force,kp,km,gammap,gammam,taup,taum,tau0,norm0,dsp,dsm,ds0);
# 	return bnd

# end

# #This is a struct which holds data defined on the boundary. First scalar
# mutable struct ScalarBndData <: AbstractBoundaryFunction
# 	U::Vector{Float64}
# 	bnd::AbstractBoundary

# 	function ScalarBndData(u::Vector{Float64},mybnd::AbstractBoundary)
# 		if ~(length(u) == mybnd.N)
# 			throw(ArgumentError("Size of input does not match boundary size"));
# 		else
# 			return new(u,mybnd)
# 		end
# 	end
# end

# #And now vector
# mutable struct VectorBndData <: AbstractBoundaryFunction
# 	U::Vector{Float64}
# 	V::Vector{Float64}
# 	bnd::AbstractBoundary

# 	function VectorBndData(u::Vector{Float64},v::Vector{Float64},mybnd::AbstractBoundary)
# 		if ~(length(u) == mybnd.N)
# 			throw(ArgumentError("Size of first input does not match boundary size"));
# 		elseif ~(length(v) == mybnd.N)
# 			throw(ArgumentError("Size of second input does not match boundary size"));
# 		else
# 			return new(u,v,mybnd)
# 		end
# 	end
# end




# ############################################
# #  Now for some spread and inter operators #
# ############################################

# function ScalarBndSpread(data::Vector{Float64},mybnd::AbstractBoundary,mygrid::AbstractGrid)
# 	if ~(length(data) == mybnd.N)
# 		throw(ArgumentError("Size of data does not match boundry"))
# 	end

# 	M = mybnd.N;
# 	Nx = mygrid.Nx;
# 	Ny = mygrid.Ny;

# 	output = zeros(Nx,Ny);

# 	for i = 1:M
# 		left = Int(floor(mybnd.bnd_loc[i,1]/mygrid.dx + 0.5));
# 		xindeces = left-1:left+2;
# 		modxin = mod1.(xindeces,Nx);
# 		coords = mygrid.dx*(xindeces .- 0.5);
# 		distance = coords .- mybnd.bnd_loc[i,1];

# 		xweights = PeskinDelta(distance/mygrid.dx)./mygrid.dx;


# 		bottom = Int(floor(mybnd.bnd_loc[i,2]/mygrid.dy + 0.5));
# 		yindeces = bottom-1:bottom+2;
# 		modyin = mod1.(yindeces,Ny);
# 		coords = mygrid.dy*(yindeces .- 0.5);
# 		distance = coords .- mybnd.bnd_loc[i,2];

# 		yweights = PeskinDelta(distance/mygrid.dy)./mygrid.dy;

# 		for (j,indx) in pairs(modxin)
# 			for (k,indy) in pairs(modyin)
# 				output[indx,indy] = output[indx,indy] + data[i]*xweights[j]*yweights[k]*mybnd.ds0[i]
# 			end
# 		end

# 	end
# 	return output
# end


# function ScalarBndInterp(data::Matrix{Float64},mybnd::AbstractBoundary,mygrid::AbstractGrid)
# 	if ~(size(data) == (mygrid.Nx,mygrid.Ny))
# 		throw(ArgumentError("Size of data does not match grid"))
# 	end

# 	M = mybnd.N;
# 	Nx = mygrid.Nx;
# 	Ny = mygrid.Ny;

# 	output = zeros(M);

# 	for i = 1:M
# 		left = Int(floor(mybnd.bnd_loc[i,1]/mygrid.dx + 0.5));
# 		xindeces = left-1:left+2;
# 		modxin = mod1.(xindeces,Nx);
# 		coords = mygrid.dx*(xindeces .- 0.5);
# 		distance = coords .- mybnd.bnd_loc[i,1];

# 		xweights = PeskinDelta(distance/mygrid.dx)./mygrid.dx;


# 		bottom = Int(floor(mybnd.bnd_loc[i,2]/mygrid.dy + 0.5));
# 		yindeces = bottom-1:bottom+2;
# 		modyin = mod1.(yindeces,Ny);
# 		coords = mygrid.dy*(yindeces .- 0.5);
# 		distance = coords .- mybnd.bnd_loc[i,2];

# 		yweights = PeskinDelta(distance/mygrid.dy)./mygrid.dy;

# 		for (j,indx) in pairs(modxin)
# 			for (k,indy) in pairs(modyin)
# 				output[i] = output[i] + data[indx,indy]*xweights[j]*yweights[k]*mygrid.dx*mygrid.dy
# 			end
# 		end

# 	end
# 	return output
# end




# function PeskinDelta(X::StepRangeLen{Float64, Float64, Float64, Int64})
# 	delta = zeros(4)
# 	delta[1] = (5.0 + 2.0*X[1] - sqrt(-7. - 12.0*X[1] - 4.0*X[1]^2) )/8.0
# 	delta[2] = (3.0 + 2.0*X[2] + sqrt(1. - 4.0*X[2] - 4.0*X[2]^2) )/8.0
# 	delta[3] = (3.0 - 2.0*X[3] + sqrt(1. + 4.0*X[3] - 4.0*X[3]^2) )/8.0
# 	delta[4] = (5.0 - 2.0*X[4] - sqrt(-7. + 12.0*X[4] - 4.0*X[4]^2) )/8.0
# 	return delta
# end



# ###########################################################
# # Simple integration techniques will be necessary as well #
# ###########################################################
# #On periodic grids simple quadrature rules should have an additional
# #Order of accuracy. Should be good enough for me. 

# function  BndIntegral(data::Vector{Float64},mybnd::AbstractBoundary)
# 	if ~(length(data) == mybnd.N)
# 		throw(ArgumentError("Size of data does not match boundry"))
# 	end

# 	summand = @. data*mybnd.ds0;
# 	output = sum(summand);
# 	return output
# end