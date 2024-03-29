# Code to be loaded into IBMProject
import FFTW

abstract type AbstractDifferentialOperator end




#Now we can define differential Operators which can act on grid data
#This one is simple and does not change the rank of the data
struct SimplePeriodicDifferentialOperator <: AbstractDifferentialOperator
	grid::PeriodicEulGrid
	Eigenvalues::Matrix{ComplexF64}
end

#This one is gradient-like. Applying it will act on scalars and produce vectors
#Inverting it will act on vectors and produce scalars
struct GradientPeriodicDifferentialOperator <: AbstractDifferentialOperator
	grid::PeriodicEulGrid
	EigenvaluesU::Matrix{ComplexF64}
	EigenvaluesV::Matrix{ComplexF64}
end

#This one is divergence-like. Applying it will act on vectors and produce scalars
#Inverting it will act on scalars and produce vectors
struct DivergencePeriodicDifferentialOperator <: AbstractDifferentialOperator
	grid::PeriodicEulGrid
	EigenvaluesU::Matrix{ComplexF64}
	EigenvaluesV::Matrix{ComplexF64}
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
    ΩX = [compfreqx[i] for i=1:mygrid.Nx,j=1:mygrid.Ny]; #put them into arrays the same size as the grid
    ΩY = [compfreqy[j] for i=1:mygrid.Nx,j=1:mygrid.Ny];

	Div = DivergencePeriodicDifferentialOperator(mygrid,ΩX,ΩY);
	return Div
end

function MakePeriodicGradient(mygrid::PeriodicEulGrid)
	ωx = FFTW.fftfreq(mygrid.Nx,mygrid.Nx/mygrid.L)
	ωy = FFTW.fftfreq(mygrid.Ny,mygrid.Ny/mygrid.H)
	compfreqx = im*ωx*(2*pi); #The eigenvalues of a single derivative in x&y
    compfreqy = im*ωy*(2*pi);
    ΩX = [compfreqx[i] for i=1:mygrid.Nx,j=1:mygrid.Ny]; #put them into arrays the same size as the grid
    ΩY = [compfreqy[j] for i=1:mygrid.Nx,j=1:mygrid.Ny];

	Grad = GradientPeriodicDifferentialOperator(mygrid,ΩX,ΩY);
	return Grad
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



#####Now the divergence-like opeators#####
#####First we do vector structs (return scalar struct)#####
function ApplyDivergenceOperator(mydata::VectorGridData,myoperator::DivergencePeriodicDifferentialOperator)
	if ~(mydata.grid === myoperator.grid)
		throw(ArgumentError("Data and operator grids do not match"));
	end
	oldUhat = FFTW.fft(mydata.U);
	oldVhat = FFTW.fft(mydata.V);
	newUhat = myoperator.EigenvaluesU.*oldUhat;
	newVhat = myoperator.EigenvaluesV.*oldVhat;
	new = real(FFTW.ifft(newUhat .+ newVhat));
	result = ScalarGridData(new,myoperator.grid);
	return result
end

#Now apply to vector data stored as two arrays
function ApplyDivergenceOperator(mydataU::Matrix{Float64},mydataV::Matrix{Float64},myoperator::DivergencePeriodicDifferentialOperator)
	if ~(size(mydataU) == (myoperator.grid.Nx,myoperator.grid.Ny))
		throw(ArgumentError("Size of first input does not match grid size"));
	elseif ~(size(mydataV) == (myoperator.grid.Nx,myoperator.grid.Ny))
		throw(ArgumentError("Size of second input does not match grid size"));
	end
	oldUhat = FFTW.fft(mydataU);
	oldVhat = FFTW.fft(mydataV);
	newUhat = myoperator.EigenvaluesU.*oldUhat;
	newVhat = myoperator.EigenvaluesV.*oldVhat;
	new = real(FFTW.ifft(newUhat .+ newVhat));
	return new
end



#####Finally, the gradient-like opeators#####
#####First we do scalar structs (return vector struct)#####
function ApplyGradientOperator(mydata::ScalarGridData,myoperator::GradientPeriodicDifferentialOperator)
	if ~(mydata.grid === myoperator.grid)
		throw(ArgumentError("Data and operator grids do not match"));
	end
	oldUhat = FFTW.fft(mydata.U);
	newUhat = myoperator.EigenvaluesU.*oldUhat;
	newVhat = myoperator.EigenvaluesV.*oldUhat;
	newU = real(FFTW.ifft(newUhat));
	newV = real(FFTW.ifft(newVhat));
	result = VectorGridData(newU,newV,myoperator.grid);
	return result
end

#Now apply to vector data stored as two arrays
function ApplyGradientOperator(mydata::Matrix{Float64},myoperator::GradientPeriodicDifferentialOperator)
	if ~(size(mydata) == (myoperator.grid.Nx,myoperator.grid.Ny))
		throw(ArgumentError("Size of first input does not match grid size"));
	end
	oldhat = FFTW.fft(mydata);
	newUhat = myoperator.EigenvaluesU.*oldhat;
	newVhat = myoperator.EigenvaluesV.*oldhat;
	newU = real(FFTW.ifft(newUhat));
	newV = real(FFTW.ifft(newVhat));
	return newU, newV
end