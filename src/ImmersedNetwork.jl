import LinearAlgebra
import StaticArrays

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
end


function MeshVertex(X::Vector{T}) where T <: Real
	if ~(length(X) == 2)
		throw(ArgumentError("Location must be a 2D location"));
	end

	point::MeshVertex = MeshVertex(X,X,0,0,false,false,[],[]);
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

	#Matrices which hold information regarding reference and current configuration
	Smat::StaticArrays.SMatrix{2,2,Float64}
	Sinv::StaticArrays.SMatrix{2,2,Float64}
	Xmat::StaticArrays.SMatrix{2,2,Float64}


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

	Smat = hcat(two.ref_coords .- one.ref_coords,three.ref_coords .- one.ref_coords);
	Xmat = hcat(two.coords .- one.coords,three.coords .- one.coords);
	Sinv = LinearAlgebra.inv(Smat);

	edge1 = MeshEdge(one,two)
	edge2 = MeshEdge(two,three)
	edge3 = MeshEdge(three,one)

	face::MeshFace = MeshFace(centroid,centroid,area,area,Smat,Sinv,Xmat,false,[edge1,edge2,edge3],[one,two,three])
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
	center = (one.coords + two.coords)./2

	edge::MeshEdge = MeshEdge(center,center,false,[one,two],[])
	return edge
end


#Lets define an Lagrangian Mesh
#And give it some constructors

mutable struct LagMesh <: AbstractNetwork
	#Number of vertices in the immersed network
	Mv::Int
	#Number of triangles in the immersed network
	Mf::Int
	#Number of edges in the immersed network
	Me::Int

	#Location of immersed boundary points, and faces that connect them
	vertices::Vector{AbstractVertex}
	faces::Vector{AbstractFace}
	edges::Vector{AbstractEdge}
end


function LagMesh(points::Vector{Vector{Float64}},triangles::Vector{Vector{S where S <: Int}})
	Mv = length(points);
	Mf = length(triangles);

	vertices = Vector{MeshVertex}(undef,Mv);
	for i = 1:Mv
		vertices[i] = MeshVertex(points[i]);
	end

	faces = Vector{MeshFace}(undef,Mf);
	for i = 1:Mf
		ipoints = triangles[i];
		faces[i] = MeshFace(vertices[ipoints[1]],vertices[ipoints[2]],vertices[ipoints[3]]);
	end

	#Euler's characteristic tells us how edges there should be, given vertices & faces
	Me = Mv + Mf - 1;
	edges = MeshEdge[];

	#We run through all the faces
	for i = 1:Mf
		#Pick out their edges
		for j = 1:3
			newedge = faces[i].edges[j];
			#If the edge hasn't been seen before, add it to the array in our mesh
			if any(x -> x.center == newedge.center,edges) == false
				push!(edges,newedge);
				println("We found a new edge")
			end
		end
	end

	#We create the mesh with all the vertices, faces & edges
	mesh::LagMesh = LagMesh(Mv,Mf,Me,vertices,faces,edges);
	#And then 'connect' it to give everyone pointers to their friends
	ConnectMesh!(mesh);
	return mesh
end


#At this point the mesh is created. Faces have edges and vertices, edges have vertices,
#Faces have both edges and vertices. 
#Now I need to write the function which connects "up" the chain. 

function ConnectMesh!(mymesh::LagMesh)

	#Lets run through all of the edges and give their vertices pointers 'home'
	for i = 1:mymesh.Me
		one,two = mymesh.edges[i].vertices
		push!(one.edges,mymesh.edges[i])
		push!(two.edges,mymesh.edges[i])
	end

	#Now we run through all the faces and give their edges AND vertices pointers 'home'.' 
	for i = 1:mymesh.Mf
		for j = 1:3
			point = mymesh.faces[i].vertices[j]
			point.area += mymesh.faces[i].area/3
			point.ref_area += mymesh.faces[i].ref_area/3
			push!(point.faces,mymesh.faces[i])

			edge = mymesh.faces[i].edges[j]
			push!(edge.faces,mymesh.faces[i])
		end
	end

	#Finally, we should run through the edges to identify boundaries
	for i = 1:mymesh.Me
		edge = mymesh.edges[i]
		if length(edge.faces) == 1
			edge.boundary = true
			edge.vertices[1].boundary = true
			edge.vertices[2].boundary = true
		end
	end
end


