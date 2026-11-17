import LinearAlgebra
import StaticArrays

abstract type AbstractVertex end
abstract type AbstractEdge end
abstract type AbstractFace end
abstract type AbstractNetwork end


#These are reserved for later version of the code
abstract type AbstractVertexFunction end
abstract type AbstractEdgeFunction end
abstract type AbstractFaceFunction end
abstract type AbstractNetworkFunction end


#Lets define the struct for a mesh vertex
#And write it a couple of constructors


mutable struct MeshVertex <: AbstractVertex
	#Location of Vertex
	coords::Point2
	ref_coords::Point2

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
	#the index of this vertex in the overall mesh
	id::Int
end

function MeshVertex(X::Vector{<: Real})
	point::MeshVertex = MeshVertex(X,0)
	return point
end

function MeshVertex(X::Vector{<: Real},id::Int)
	if ~(length(X) == 2)
		throw(ArgumentError("Location must be a 2D location"));
	end

	point::MeshVertex = MeshVertex(Point2(X),Point2(X),0,0,false,false,[],[],id);
	return point
end


#Now lets define the struct for a mesh face
#And give it some construtors

mutable struct MeshFace <: AbstractFace
	#Location of face centroid
	centroid::Point2
	ref_centroid::Point2

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

function MeshFace(vert1::MeshVertex,vert2::MeshVertex,vert3::MeshVertex)
	centroid = TriangleCentroid(vert1.coords,vert2.coords,vert3.coords);
	area = TriangleArea(vert1.coords,vert2.coords,vert3.coords);

	Smat = TriangleEdgeMatrix(vert1.ref_coords,vert2.ref_coords,vert3.ref_coords);
	Xmat = TriangleEdgeMatrix(vert1.coords,vert2.coords,vert3.coords);
	Sinv = LinearAlgebra.inv(Smat);

	edge1 = MeshEdge(vert1,vert2)
	edge2 = MeshEdge(vert2,vert3)
	edge3 = MeshEdge(vert3,vert1)

	face::MeshFace = MeshFace(centroid,centroid,area,area,Smat,Sinv,Xmat,false,[edge1,edge2,edge3],[vert1,vert2,vert3])
	return face
end

function MeshFace(vertices::Vector{MeshVertex})
	if ~(length(vertices) == 3)
		throw(ArgumentError("Face must be defined by 3 vertices"));
	end
	vert1 = vertices[1];
	vert2 = vertices[2];
	vert3 = vertices[3];
	face::MeshFace = MeshFace(vert1,vert2,vert3)
	return face
end


#Now lets define the struct for a mesh edge
#And give it some construtors

mutable struct MeshEdge <: AbstractEdge
	#Location of face centroid
	center::Point2
	ref_center::Point2

	#Booleans that indicate if this is a boundary face
	boundary::Bool

	#Now we need vectors of this nodes faces and vertices
	vertices::Vector{AbstractVertex}
	faces::Vector{AbstractFace}

end

function MeshEdge(vert1::MeshVertex,vert2::MeshVertex)
	center = Midpoint(vert1.coords,vert2.coords);

	edge::MeshEdge = MeshEdge(center,center,false,[vert1,vert2],[])
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


function LagMesh(points::Vector{<:Vector{<:Real}},triangles::Vector{<:Vector{<:Integer}})
	Mv = length(points);
	Mf = length(triangles);

	vertices = Vector{MeshVertex}(undef,Mv);
	for i = 1:Mv
		vertices[i] = MeshVertex(points[i],i);
	end

	faces = Vector{MeshFace}(undef,Mf);
	for i = 1:Mf
		ipoints = triangles[i];
		faces[i] = MeshFace(vertices[ipoints[1]],vertices[ipoints[2]],vertices[ipoints[3]]);
	end

	#Create an empty dictionary of keys & edges
	edgedict = Dict{Tuple{Int,Int},MeshEdge}()
	#Create an empty array to hold the unique edges
	edges = MeshEdge[];

	#Loop through the faces
	for i = 1:Mf
	    ip = triangles[i];
	    #Loop through the points of this face
	    for j = 1:3
	    	#Grab the j'th point in this face, plus the "next" one
	        a = Int(ip[j]);
	        b = Int(ip[mod1(j+1,3)]);
	        #Turn those integers into an _ordered_ tuple
	        key = (min(a,b), max(a,b));
	        #Check if this tuple is in the dictionary yet. 
	        if haskey(edgedict, key)
	            # Seen before: point this face at the shared edge object
	            faces[i].edges[j] = edgedict[key];
	        else
	            # First time: this face's edge becomes the canonical one
	            edgedict[key] = faces[i].edges[j];
	            push!(edges, faces[i].edges[j]);
	        end
	    end
	end

	#Chech the final negth of the list of edges
	Me = length(edges); 

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
		vert1,vert2 = mymesh.edges[i].vertices
		push!(vert1.edges,mymesh.edges[i])
		push!(vert2.edges,mymesh.edges[i])
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
			edge.faces[1].boundary = true
		end
	end
end


