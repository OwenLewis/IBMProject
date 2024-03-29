# Code to be loaded into IBMProject

abstract type AbstractGrid end
abstract type AbstractGridFunction end

#Lets define an Eulerian Grid
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

#If data is specified, check size
function ScalarGridData(ux::Matrix{T},mygrid::AbstractGrid) where T<:Real
	if ~(size(ux) == size(mygrid.X))
		throw(ArgumentError(ux,"Size of input does not match grid size"));
	else
		griddata::ScalarGridData = ScalarGridData(ux,mygrid);
		return griddata
	end
end

#If no data is specified, initialize with zeros
function VectorGridData(mygrid::AbstractGrid)
	U = zeros(mygrid.Nx,mygrid.Ny);
	V = zeros(mygrid.Nx,mygrid.NY)
	griddata::ScalarGridData = VectorGridData(U,V,mygrid);
	return griddata
end

#If data is specified, check size
function VectorGridData(ux::Matrix{T},uy::Matrix{T},mygrid::AbstractGrid) where T<:Real
	if ~(size(ux) == size(mygrid.X))
		throw(ArgumentError(ux,"Size of input does not match grid size"));
	elseif ~(size(uy) == size(mygrid.X))
		throw(ArgumentError("Size of second input does not match grid size"));
	else
		griddata::ScalarGridData = VectorGridData(ux,uy,mygrid);
		return griddata
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