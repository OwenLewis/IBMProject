# Code to be loaded into IBMProject
import FFTW

abstract type AbstractDifferentialOperator end




#Now we can define differential Operators which can act on grid data
#This one is simple and does not change the rank of the data
struct SimplePeriodicDifferentialOperator <: AbstractDifferentialOperator
	grid::PeriodicEulGrid
	Eigenvalues::Matrix{Float64}
end

#This one is gradient-like. Applying it will act on scalars and produce vectors
#Inverting it will act on vectors and produce scalars
struct GradientPeriodicDifferentialOperator <: AbstractDifferentialOperator
	grid::PeriodicEulGrid
	EigenvaluesU::Matrix{Float64}
	EigenvaluesV::Matrix{Float64}
end

#This one is divergence-like. Applying it will act on vectors and produce scalars
#Inverting it will act on scalars and produce vectors
struct DivergencePeriodicDifferentialOperator <: AbstractDifferentialOperator
	grid::PeriodicEulGrid
	EigenvaluesU::Matrix{Float64}
	EigenvaluesV::Matrix{Float64}
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
	Lap = SimplePeriodicDifferentialOperator(mygrid,applyeigs);
	return Lap
end

function MakePeriodicDivergence(mygrid::PeriodicEulGrid)
	ωx = FFTW.fftfreq(mygrid.Nx,mygrid.Nx/mygrid.L)
	ωy = FFTW.fftfreq(mygrid.Ny,mygrid.Ny/mygrid.H)
	compfreqx = im*ωx*(2*pi); #The eigenvalues of a single derivative in x&y
    compfreqy = im*ωy*(2*pi);
    ΩX = (compfreqx[i] for i=1:mygrid.Nx,j=1:mygrid.Ny); #put them into arrays the same size as the grid
    ΩY = (compfreqy[j] for i=1:mygrid.Nx,j=1:mygrid.Ny);

    applyeigsU = ΩX;
    applyeigsV = ΩY;

	Div = DivergencePeriodicDifferentialOperator(mygrid,applyeigsU,applyeigsV);
	return Div
end



###############################################
#   This is how we apply/invert operators     #
###############################################


#####Simple ones first, there are a few cases#####
#####First we do scalar data, either struct or simple arrays#####
#This applys a simple differential operator to a scalar function on a periodic grid
function ApplySimpleOperator(mydata::ScalarGridData,myoperator::SimplePeriodicDifferentialOperator)
	if ~(mydata.grid === myoperator.grid)
		throw(ArgumentError("Data and operator grids do not match"));
	end
	oldhat = FFTW.fft(mydata.U);
	newhat = myoperator.Eigenvalues.*oldhat;
	new = real(FFTW.ifft(newhat));
	result = ScalarGridData(new,myoperator.grid);
	return result
end

#Same, but no need for a struct. Will apply to an array of data
function ApplySimpleOperator(mydata::Matrix{Float64},myoperator::SimplePeriodicDifferentialOperator)
	if ~(size(mydata) == (myoperator.grid.Nx,myoperator.grid.Ny))
		throw(ArgumentError("Size of input does not match grid size of operator"));
	end
	oldhat = FFTW.fft(mydata);
	newhat = myoperator.Eigenvalues.*oldhat;
	new = real(FFTW.ifft(newhat));
	return new
end


####Now the same as above, but for vector structs or pairs of arrays####

#Applying the operator to a vector grid data struct
function ApplySimpleOperator(mydata::VectorGridData,myoperator::SimplePeriodicDifferentialOperator)
	if ~(mydata.grid === myoperator.grid)
		throw(ArgumentError("Data and operator grids do not match"));
	end
	oldUhat = FFTW.fft(mydata.U);
	oldVhat = FFTW.fft(mydata.V);
	newUhat = myoperator.Eigenvalues.*oldUhat;
	newVhat = myoperator.Eigenvalues.*oldVhat;
	newU = real(FFTW.ifft(newUhat));
	newV = real(FFTW.ifft(newVhat));
	result = VectorGridData(newU,newV,myoperator.grid);
	return result
end

#Applying it to two arrays
function ApplySimpleOperator(mydataU::Matrix{Float64},mydataV::Matrix{Float64},myoperator::SimplePeriodicDifferentialOperator)
	if ~(size(mydataU) == (myoperator.grid.Nx,myoperator.grid.Ny))
		throw(ArgumentError("Size of first input does not match grid size"));
	elseif ~(size(mydataV) == (myoperator.grid.Nx,myoperator.grid.Ny))
		throw(ArgumentError("Size of second input does not match grid size"));
	end
	oldUhat = FFTW.fft(mydataU);
	oldVhat = FFTW.fft(mydataV);
	newUhat = myoperator.Eigenvalues.*oldUhat;
	newVhat = myoperator.Eigenvalues.*oldVhat;
	newU = real(FFTW.ifft(newUhat));
	newV = real(FFTW.ifft(newVhat));
	return newU,newV
end
