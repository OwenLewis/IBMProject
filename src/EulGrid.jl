# Code related to the construction of Eulerian grids
# and data types to store functions on said grids.

import StaticArrays 
import LinearAlgebra

abstract type AbstractGrid end
abstract type AbstractGridFunction end

#Lets define an Eulerian Grid
#NOTE ON GRID IDENTITY: the operators in DiffOpers.jl check that data and operator
#share a grid using `===`. For this struct that compares the numeric fields by value
#but X and Y by object identity, so two grids built by separate calls to
#PeriodicEulGrid(...) with identical arguments are NOT `===` (and not `==` either,
#since `==` falls back to `===`). Build one grid and reuse that object everywhere.
#If that ever becomes a nuisance, define a SameGrid(g1,g2) helper comparing
#Lx,Ly,Ux,Uy,Nx,Ny (the rest is derived) and use it in the Apply*Operator checks.
struct PeriodicEulGrid <: AbstractGrid
	Lx::Float64
	Ly::Float64
	Ux::Float64
	Uy::Float64
	Nx::Int
	Ny::Int
	dx::Float64
	dy::Float64
	X::Matrix{Float64}
	Y::Matrix{Float64}
	L::Float64
	H::Float64
end

#Lets define scalar data on an Eulerian Grid (with a constructor for exception checking)
mutable struct ScalarGridData <: AbstractGridFunction
	U::Matrix{Float64}
	grid::AbstractGrid

	function ScalarGridData(ux::Matrix{Float64},mygrid::AbstractGrid)
		if ~(size(ux) == size(mygrid.X))
			throw(ArgumentError("Size of input does not match grid size"));
		else
			return new(ux,mygrid)
		end
	end
end

#Lets define vector data on the Eulerian Grid
mutable struct VectorGridData <: AbstractGridFunction
	U::Matrix{Float64}
	V::Matrix{Float64}
	grid::AbstractGrid

	function VectorGridData(ux::Matrix{Float64},uy::Matrix{Float64},mygrid::AbstractGrid)
		if ~(size(ux) == size(mygrid.X))
			throw(ArgumentError("Size of first input does not match grid size"));
		elseif ~(size(uy) == size(mygrid.X))
			throw(ArgumentError("Size of second input does not match grid size"));
		else
			return new(ux,uy,mygrid)
		end 
	end
end

###############################################
#          Now some grid constructors         #
###############################################

#No need to specify the redundant information
function PeriodicEulGrid(Lx::T,Ly::T,Ux::T,Uy::T,Nx::Int,Ny::Int) where T<:Real
	L = Ux - Lx;
	H = Uy - Ly;
	dx = L/Nx;
	dy = H/Ny;
	X = [(i-1/2)*dx + Lx for i=1:Nx, j = 1:Ny];
	Y = [(j-1/2)*dy + Ly for i=1:Nx, j = 1:Ny];
	grid::PeriodicEulGrid = PeriodicEulGrid(Lx,Ly,Ux,Uy,Nx,Ny,dx,dy,X,Y,L,H);
	return grid
end

#If lower bounds are not specified, we assume they're zero
function PeriodicEulGrid(Ux::T,Uy::T,Nx::Int,Ny::Int) where T<:Real
	Lx = 0.0;
	Ly = 0.0;
	L = Ux - Lx;
	H = Uy - Ly;
	dx = L/Nx;
	dy = H/Ny;
	X = [(i-1/2)*dx + Lx for i=1:Nx, j = 1:Ny];
	Y = [(j-1/2)*dy + Ly for i=1:Nx, j = 1:Ny];
	grid::PeriodicEulGrid = PeriodicEulGrid(Lx,Ly,Ux,Uy,Nx,Ny,dx,dy,X,Y,L,H);
	return grid
end

#If upper bounds are not specified, we assume they're unity
function PeriodicEulGrid(Nx::Int,Ny::Int)
	Lx = 0.0;
	Ly = 0.0;
	Ux = 1.0;
	Uy = 1.0;
	L = Ux - Lx;
	H = Uy - Ly;
	dx = L/Nx;
	dy = H/Ny;
	X = [(i-1/2)*dx + Lx for i=1:Nx, j = 1:Ny];
	Y = [(j-1/2)*dy + Ly for i=1:Nx, j = 1:Ny];
	grid::PeriodicEulGrid = PeriodicEulGrid(Lx,Ly,Ux,Uy,Nx,Ny,dx,dy,X,Y,L,H);
	return grid
end

###############################################
#        Now some grid data constructors      #
###############################################
#If no data is specified, initialize with zeros
function ScalarGridData(mygrid::AbstractGrid)
	U = zeros(mygrid.Nx,mygrid.Ny);
	griddata::ScalarGridData = ScalarGridData(U,mygrid);
	return griddata
end

#If data is specified, convert type
function ScalarGridData(ux::Matrix{T},mygrid::AbstractGrid) where T<:Real
	griddata::ScalarGridData = ScalarGridData(convert(Matrix{Float64}, ux),mygrid)
	return griddata
end

#If no data is specified, initialize with zeros
function VectorGridData(mygrid::AbstractGrid)
	U = zeros(mygrid.Nx,mygrid.Ny);
	V = zeros(mygrid.Nx,mygrid.Ny)
	griddata::VectorGridData = VectorGridData(U,V,mygrid);
	return griddata
end

#If data is specified, convert type
function VectorGridData(ux::Matrix{T},uy::Matrix{S},mygrid::AbstractGrid) where {T<:Real, S<:Real}
	griddata::VectorGridData = VectorGridData(convert(Matrix{Float64}, ux),convert(Matrix{Float64}, uy),mygrid);
	return griddata
end

# Create a vector field using two Scalar fields
function VectorGridData(a::ScalarGridData,b::ScalarGridData)
	CheckSameGrid(a,b,"concatenated")
	return VectorGridData(copy(a.U),copy(b.U),a.grid);
end



#################################################################
#  Wrappers to accomplish Matrix{Point2} <---> VectorGridData   #
#################################################################

#Get a vector field from a matrix of Point2
function VectorGridData(M::Matrix{<:StaticArrays.SVector{2,<:Real}},mygrid::AbstractGrid)
	return VectorGridData(first.(M),last.(M),mygrid)
end

# Get a matrix of Point2 from a vector field
function Point2Matrix(v::VectorGridData) 
	return Point2.(v.U,v.V);
end

# Some helpful methods to get/assign grid data values at various indices
# Not sure if these will ever be helpful, but they aren't likely to hurt

# Return v[i,j] at one cell as a Point2
Base.getindex(v::VectorGridData,i::Integer,j::Integer) = Point2(v.U[i,j],v.V[i,j])

# Assign v[i,j] = p using a 2D vector at one cell
function Base.setindex!(v::VectorGridData,p::StaticArrays.SVector{2,<:Real},i::Integer,j::Integer)
	v.U[i,j] = p[1];
	v.V[i,j] = p[2];
	return v
end

# Accumulate into one cell; this is what spreading needs 
#(a plain v[i,j] += f also works through the two methods above, but allocates less this way)
function AddAt!(v::VectorGridData,i::Integer,j::Integer,f::StaticArrays.SVector{2,<:Real})
	v.U[i,j] += f[1];
	v.V[i,j] += f[2];
	return v
end

# Helper function to return one component of a vector field as a scalar field
function Component(v::VectorGridData,k::Integer)
	if k == 1
		return ScalarGridData(copy(v.U),v.grid);
	elseif k == 2
		return ScalarGridData(copy(v.V),v.grid);
	else
		throw(ArgumentError("Vector grid data has components 1 and 2, got $k."));
	end
end

###############################################
#      Lets define some grid data algebra     #
###############################################

# Adding & subtracting scalar fields
function Base.:+(a::ScalarGridData,b::ScalarGridData)
	CheckSameGrid(a,b,"added")
	return ScalarGridData(a.U .+ b.U,a.grid)
end
function Base.:-(a::ScalarGridData,b::ScalarGridData)
	CheckSameGrid(a,b,"subtracted")
	return ScalarGridData(a.U .- b.U,a.grid)
end

#Adding and subtracting vector fields
function Base.:+(a::VectorGridData,b::VectorGridData)
	CheckSameGrid(a,b,"added")
	return VectorGridData(a.U .+ b.U,a.V .+ b.V,a.grid)
end
function Base.:-(a::VectorGridData,b::VectorGridData)
	CheckSameGrid(a,b,"subtracted")
	return VectorGridData(a.U .- b.U,a.V .- b.V,a.grid)
end

# Scalar operations on scalar fields
function Base.:-(a::ScalarGridData) 
	return ScalarGridData(-a.U,a.grid);
end
function Base.:*(c::Real,a::ScalarGridData) 
	return ScalarGridData(c.*a.U,a.grid);
end
function Base.:*(a::ScalarGridData,c::Real)
	return c*a;
end
function Base.:/(a::ScalarGridData,c::Real) 
	return ScalarGridData(a.U./c,a.grid);
end

#Scalar operations on vector fields
function Base.:-(a::VectorGridData)
	return VectorGridData(-a.U,-a.V,a.grid);
end
function  Base.:*(c::Real,a::VectorGridData) 
	return VectorGridData(c.*a.U,c.*a.V,a.grid);
end
function Base.:*(a::VectorGridData,c::Real) 
	return c*a;
end
function  Base.:/(a::VectorGridData,c::Real)
	return VectorGridData(a.U./c,a.V./c,a.grid);
end

#Pointwise mult of scalar/vector fields
function Base.:*(a::ScalarGridData,b::ScalarGridData)
	CheckSameGrid(a,b,"multiplied")
	return ScalarGridData(a.U.*b.U,a.grid);
end
function Base.:/(a::ScalarGridData,b::ScalarGridData)
	CheckSameGrid(a,b,"divided")
	return ScalarGridData(a.U./b.U,a.grid);
end
function Base.:*(a::ScalarGridData,b::VectorGridData) 
	CheckSameGrid(a,b,"multiplied")
	return VectorGridData(a.U.*b.U,a.U.*b.V,a.grid);
end
function Base.:*(b::VectorGridData,a::ScalarGridData)
	CheckSameGrid(a,b,"multiplied")
	return a*b;
end

#Pointwise products of two vector fields
#This is a dot-product. It returns a scalar field.
function LinearAlgebra.dot(a::VectorGridData,b::VectorGridData)
	CheckSameGrid(a,b,"dotted")
	return ScalarGridData(a.U.*b.U + a.V.*b.V,a.grid);
end
#This is a 2D cross product. It currently returns a scalar field
#Things will need to change substantially if we move to 3D
function LinearAlgebra.cross(a::VectorGridData,b::VectorGridData)
	CheckSameGrid(a,b,"crossed")
	return ScalarGridData(a.U.*b.V - a.V.*b.U,a.grid);
end



# Helper Function to return a boolean determining if two grid functions are on the same grid
function CheckSameGrid(A::ScalarGridData,B::ScalarGridData,what::String)
	if ~(A.grid === B.grid)
		throw(ArgumentError("Scalar fields to be $what are on different grids"));
	end
end
function CheckSameGrid(A::ScalarGridData,B::VectorGridData,what::String)
	if ~(A.grid === B.grid)
		throw(ArgumentError("Scalar field and vector field to be $what are on different grids"));
	end
end
function CheckSameGrid(A::VectorGridData,B::ScalarGridData,what::String)
	if ~(A.grid === B.grid)
		throw(ArgumentError("Vector field and scalar field to be $what are on different grids"));
	end
end
function CheckSameGrid(A::VectorGridData,B::VectorGridData,what::String)
	if ~(A.grid === B.grid)
		throw(ArgumentError("Vector fields to be $what are on different grids"));
	end
end



###########################################################
# Simple integration techniques will be necessary as well #
###########################################################
#On periodic grids simple quadrature rules should have an additional
#Order of accuracy. Should be good enough for me. 

function  GridIntegral(data::Matrix{Float64},mygrid::AbstractGrid)
	if ~(size(data) == (mygrid.Nx,mygrid.Ny))
		throw(ArgumentError("Size of data does not match grid"))
	end

	summand = @. data*mygrid.dx*mygrid.dy;
	output = sum(summand);
	return output
end