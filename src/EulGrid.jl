# Code to be loaded into IBMProject
import FFTW

abstract type AbstractGrid end
abstract type AbstractGridFunction end
abstract type AbstractDifferentialOperator end

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

#Now we can define differential Operators which can act on grid data
struct PeriodicDifferentialOperator <: AbstractDifferentialOperator
	grid::PeriodicEulGrid
	applyEigenvalues::Matrix{Float64}
	invertEigenvalues::Matrix{Float64}
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

###############################################
#         Now for some useful functions       #
###############################################


#This is to make a periodic laplacian of type PeriodicDifferentialOperator
function MakePeriodicLaplacian(mygrid::PeriodicEulGrid)
	ωx = FFTW.fftfreq(mygrid.Nx,mygrid.Nx/mygrid.L)
	ωy = FFTW.fftfreq(mygrid.Ny,mygrid.Ny/mygrid.H)
	compfreqx = im*ωx*(2*pi); #The eigenvalues of a single derivative in x&y
    compfreqy = im*ωy*(2*pi);
    ΩX = (compfreqx[i] for i=1:mygrid.Nx,j=1:mygrid.Ny); #put them into arrays the same size as the grid
    ΩY = (compfreqy[j] for i=1:mygrid.Nx,j=1:mygrid.Ny);

    applyeigs = (ΩX.^2 + ΩY.^2); #Eigenvalues of the laplacian
    inteigs = 1 ./applyeigs; #The inverse of the eigenvalues of the laplacian
    inteigs[1,1] = 0;    #Zero out the 0-0 eigenvalue
	Lap = PeriodicDifferentialOperator(mygrid,applyeigs,inteigs);
	return Lap
end


#This function applys a differential operator to a scalar function on a periodic grid
function ApplySingleOperator(mydata::ScalarGridData,myoperator::PeriodicDifferentialOperator,mygrid::PeriodicEulGrid)
	oldhat = FFTW.fft(mydata.U);
	newhat = myoperator.applyEigenvalues.*oldhat;
	new = real(FFTW.ifft(newhat));
	result = ScalarGridData(new,mygrid);
	return result
end

#This function inverts a differential operator on a scalar function on a grid
function InvertSingleOperator(mydata::ScalarGridData,myoperator::PeriodicDifferentialOperator,mygrid::PeriodicEulGrid)
	oldhat = FFTW.fft(mydata.U);
	newhat = myoperator.invertEigenvalues.*oldhat;
	new = real(FFTW.ifft(newhat));
	result = ScalarGridData(new,mygrid);
	return result
end

#Now we need two new methods for the above function, but for vector data. 
function ApplySingleOperator(mydata::VectorGridData,myoperator::PeriodicDifferentialOperator,mygrid::PeriodicEulGrid)
	oldUhat = FFTW.fft(mydata.U);
	oldVhat = FFTW.fft(mydata.V);
	newUhat = myoperator.applyEigenvalues.*oldUhat;
	newVhat = myoperator.applyEigenvalues.*oldVhat;
	newU = real(FFTW.ifft(newUhat));
	newV = real(FFTW.ifft(newVhat));
	result = ScalarGridData(newU,newV,mygrid);
	return result
end

function InvertSingleOperator(mydata::VectorGridData,myoperator::PeriodicDifferentialOperator,mygrid::PeriodicEulGrid)
	oldUhat = FFTW.fft(mydata.U);
	oldVhat = FFTW.fft(mydata.V);
	newUhat = myoperator.invertEigenvalues.*oldUhat;
	newVhat = myoperator.invertEigenvalues.*oldVhat;
	newU = real(FFTW.ifft(newUhat));
	newV = real(FFTW.ifft(newVhat));
	result = ScalarGridData(newU,newV,mygrid);
	return result
end